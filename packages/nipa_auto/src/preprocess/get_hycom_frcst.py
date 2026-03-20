
#!/usr/bin/env python3
from __future__ import annotations

import os
import time
import re
import yaml
import numpy as np
import xarray as xr
from datetime import datetime as dt, timedelta
from typing import Any, Dict, Tuple, Optional
from pathlib import Path

BASE = Path(__file__).resolve().parent
cfg_path = BASE / "config.yaml"

# =========================
# USER TUNABLE
# =========================
FMRC_URL = "https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_all/FMRC_ESPC-D-V02_all_best.ncd"

# 저장 샘플링 간격(시간)
SAVE_STEP_HOURS = 3

# ROMS-style time reference
TIME_REF = dt(2000, 1, 1, 0, 0, 0)

# nearest 허용 오차(분)
NEAREST_TOL_MINUTES = 180


# =========================
# helpers
# =========================
def as_dt(x: Any) -> dt:
    if isinstance(x, dt):
        return x
    if isinstance(x, str):
        return dt.strptime(x, "%Y-%m-%d %H:%M:%S")
    raise TypeError(f"Unsupported time type: {type(x)} -> {x}")

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def pick_run_date(cfg: Dict[str, Any], t_st: dt, t_ed: dt) -> str:
    cand_keys = ["initdate", "init_date", "run_date", "RUN_DATE", "RUN_DATE_ISO"]
    for k in cand_keys:
        if k in cfg and isinstance(cfg[k], str):
            s = cfg[k].strip()
            for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
                try:
                    return dt.strptime(s, fmt).strftime("%Y%m%d")
                except Exception:
                    pass
    mid = t_st + (t_ed - t_st) / 2
    return mid.strftime("%Y%m%d")

def detect_latlon_names(ds: xr.Dataset) -> Tuple[str, str]:
    for cand in ("lat", "latitude", "Latitude", "LAT"):
        if cand in ds.coords or cand in ds.variables:
            lat_name = cand
            break
    else:
        raise KeyError("Cannot find latitude coordinate (lat/latitude/Latitude/LAT).")

    for cand in ("lon", "longitude", "Longitude", "LON"):
        if cand in ds.coords or cand in ds.variables:
            lon_name = cand
            break
    else:
        raise KeyError("Cannot find longitude coordinate (lon/longitude/Longitude/LON).")

    return lat_name, lon_name

def safe_slice(ds: xr.Dataset, lat_name: str, lon_name: str,
               lat_rng: Tuple[float, float], lon_rng: Tuple[float, float]) -> xr.Dataset:
    lat0, lat1 = lat_rng
    lon0, lon1 = lon_rng

    if ds[lat_name].size >= 2 and float(ds[lat_name][0]) > float(ds[lat_name][-1]):
        ds = ds.sortby(lat_name)
    if ds[lon_name].size >= 2 and float(ds[lon_name][0]) > float(ds[lon_name][-1]):
        ds = ds.sortby(lon_name)

    return ds.sel({lat_name: slice(lat0, lat1), lon_name: slice(lon0, lon1)})

def make_regular_times(t_st: dt, t_ed: dt, step_hours: int) -> np.ndarray:
    step = timedelta(hours=step_hours)
    out = []
    cur = t_st
    while cur <= t_ed:
        out.append(np.datetime64(cur))
        cur += step
    return np.array(out)

def open_with_retry(url: str, tries: int = 5, wait_sec: int = 30) -> xr.Dataset:
    last = None
    for i in range(1, tries + 1):
        try:
            return xr.open_dataset(url, decode_times=False)
        except Exception as e:
            last = e
            if i == tries:
                break
            print(f"[WARN] open_dataset failed ({i}/{tries}): {e} -> sleep {wait_sec}s", flush=True)
            time.sleep(wait_sec)
    raise RuntimeError(f"Failed to open dataset after {tries} tries: {last}")

def _parse_ref_datetime(ref: str) -> dt:
    s = ref.strip()
    toks = s.split()
    if toks and toks[-1].upper() in ("UTC", "GMT", "Z"):
        toks = toks[:-1]
    s = " ".join(toks)

    if " " in s:
        dpart, tpart = s.split(" ", 1)
        if "." in tpart:
            tpart = tpart.split(".", 1)[0]
        s = f"{dpart} {tpart}".strip()

    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return dt.strptime(s, fmt)
        except Exception:
            pass
    raise ValueError(f"Cannot parse ref datetime: {ref!r} -> {s!r}")

def decode_time_coord(ds: xr.Dataset, tname: str) -> xr.Dataset:
    tv = ds[tname]
    units = tv.attrs.get("units", "")
    if "since" not in units:
        raise ValueError(f"Unexpected units for {tname}: {units!r}")

    unit, _, ref = units.partition(" since ")
    unit = unit.strip().lower()
    ref_dt = _parse_ref_datetime(ref)

    vals = np.asarray(tv.values, dtype=float)

    if unit.startswith("hour"):
        times = np.array([np.datetime64(ref_dt + timedelta(hours=float(v))) for v in vals])
    elif unit.startswith("day"):
        times = np.array([np.datetime64(ref_dt + timedelta(days=float(v))) for v in vals])
    elif unit.startswith("min"):
        times = np.array([np.datetime64(ref_dt + timedelta(minutes=float(v))) for v in vals])
    elif unit.startswith("sec"):
        times = np.array([np.datetime64(ref_dt + timedelta(seconds=float(v))) for v in vals])
    else:
        raise ValueError(f"Unsupported unit for {tname}: {units!r}")

    return ds.assign_coords({tname: times})

def decode_all_time_coords(ds: xr.Dataset) -> xr.Dataset:
    """
    time, time1, time2... 중 units에 'since'가 들어있는 coord만 datetime64로 디코드.
    (hours since analysis 같은 비표준은 건너뜀)
    """
    for cname in list(ds.coords):
        if not cname.lower().startswith("time"):
            continue
        units = ds[cname].attrs.get("units", "")
        if "since" in units:
            try:
                ds = decode_time_coord(ds, cname)
            except Exception:
                pass
    return ds

def subset_vars(ds: xr.Dataset) -> xr.Dataset:
    want = [v for v in ("surf_el", "water_temp", "salinity", "water_u", "water_v") if v in ds.data_vars]
    if not want:
        raise KeyError("None of required HYCOM vars found: surf_el, water_temp, salinity, water_u, water_v")
    return ds[want]

def find_time_dim(da: xr.DataArray) -> Optional[str]:
    for d in da.dims:
        if d.lower().startswith("time"):
            return d
    return None

def resample_var_to_target(
    da: xr.DataArray,
    t_st: dt, t_ed: dt,
    target_times: np.ndarray,
    tol: np.timedelta64
) -> xr.DataArray:
    """
    변수마다 time dim이 time/time1/...로 다를 수 있으니:
      - 그 time dim 기준 window slice
      - target_times로 nearest 선택
      - dim 이름을 'time'으로 통일
    """
    tdim = find_time_dim(da)
    if tdim is None:
        return da
    if tdim not in da.coords:
        raise ValueError(f"{da.name}: has time dim {tdim} but no coord")

    da2 = da.sel({tdim: slice(np.datetime64(t_st), np.datetime64(t_ed))})
    da2 = da2.sel({tdim: target_times}, method="nearest", tolerance=tol)

    if tdim != "time":
        da2 = da2.rename({tdim: "time"})
    return da2


# ===== units normalize for saved file (optional safety) =====
def normalize_time_units(units: str) -> str:
    """
    'hours since 2026-01-13 12:00:00.000 UTC' -> 'hours since 2026-01-13 12:00:00'
    """
    if not isinstance(units, str):
        return units
    s = units.strip()
    s = re.sub(r"\s+(UTC|GMT|Z)\s*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"(\d{2}:\d{2}:\d{2})\.\d+", r"\1", s)
    return s

def normalize_all_time_coord_units(ds: xr.Dataset) -> xr.Dataset:
    for cname in list(ds.coords):
        if cname.lower().startswith("time"):
            u = ds[cname].attrs.get("units", None)
            if isinstance(u, str) and "since" in u:
                nu = normalize_time_units(u)
                ds[cname].attrs["units"] = nu
                ds[cname].encoding["units"] = nu
    return ds


# =========================
# main
# =========================
def main() -> int:
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t_st = as_dt(cfg["bry_start_date"])
    t_ed = as_dt(cfg["bry_end_date"])
    lat_rng = tuple(cfg["region"]["lat"])
    lon_rng = tuple(cfg["region"]["lon"])

    base_dir = cfg["output"]["base_dir"]
    out_dir = os.path.join(base_dir, "hycom")
    ensure_dir(out_dir)

    run_yyyymmdd = pick_run_date(cfg, t_st, t_ed)
    out_path = os.path.join(out_dir, f"hycom_korea_{run_yyyymmdd}.nc")

    print("[INFO] HYCOM FMRC url:", FMRC_URL)
    print("[INFO] window:", t_st, "to", t_ed)
    print("[INFO] save step (hours):", SAVE_STEP_HOURS)
    print("[INFO] region lat:", lat_rng, "lon:", lon_rng)
    print("[INFO] output:", out_path)

    print("[INFO] opening dataset (decode_times=False)...")
    ds = open_with_retry(FMRC_URL, tries=5, wait_sec=30)

    lat_name, lon_name = detect_latlon_names(ds)
    print("[INFO] lat/lon coords:", {"lat": lat_name, "lon": lon_name})

    print("[INFO] decoding all time-like coords (time/time1/...) ...")
    ds = decode_all_time_coords(ds)

    ds = subset_vars(ds)

    # 공간 subset
    ds = safe_slice(ds, lat_name, lon_name, lat_rng, lon_rng)

    # target times
    target_times = make_regular_times(t_st, t_ed, SAVE_STEP_HOURS)
    tol = np.timedelta64(int(NEAREST_TOL_MINUTES), "m")

    # 변수별 time dim으로 resample 후 time 통일
    print("[INFO] resampling each variable to target_times and unifying time dim ...")
    out_vars = {}
    for v in ds.data_vars:
        da = ds[v]
        da2 = resample_var_to_target(da, t_st, t_ed, target_times, tol)
        out_vars[v] = da2

        tdim = find_time_dim(da)
        if tdim is not None:
            print(f"  - {v}: {tdim} -> time, nt={da2.sizes.get('time', 'NA')}")

    out = xr.Dataset(out_vars)

    # (안전) time 말고 time1/time2 남았으면 제거
    drop_coords = [c for c in out.coords if c.lower().startswith("time") and c != "time"]
    if drop_coords:
        out = out.drop_vars(drop_coords)

    # time을 ROMS 스타일 days since로 변환
    print("[INFO] converting unified time coord to days since 2000-01-01 ...")
    tvals = out["time"].values.astype("datetime64[s]")
    secs = (tvals - np.datetime64(TIME_REF)).astype("timedelta64[s]").astype(np.int64)
    days = secs / 86400.0
    out = out.assign_coords(time=("time", days))
    out["time"].attrs["units"] = "days since 2000-01-01 00:00:00"
    out["time"].attrs["calendar"] = "proleptic_gregorian"

    # 저장 시 units 문자열 정규화(UTC/.000 제거)
    out = normalize_all_time_coord_units(out)
    out["time"].encoding["units"] = out["time"].attrs["units"]

    enc = {v: {"zlib": True, "complevel": 4} for v in out.data_vars}

    tmp = out_path + ".tmp"
    print("[INFO] writing netcdf ...")
    out.to_netcdf(tmp, encoding=enc, unlimited_dims=["time"])
    os.replace(tmp, out_path)

    print("[OK] saved:", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
