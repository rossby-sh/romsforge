#!/usr/bin/env python3
from __future__ import annotations

import os
import re
from pathlib import Path
from datetime import datetime as dt, timedelta
from typing import Any, Dict, List, Tuple, Optional

import requests
import subprocess
import numpy as np
import xarray as xr
import yaml


BASE = Path(__file__).resolve().parent
cfg_path = BASE / "config.yaml"

S3_BASE = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
RES = "0p25"
RUN_HOURS = (0, 6, 12, 18)

STEP_HOURS = 3
TIME_REF = dt(2000, 1, 1, 0, 0, 0)

HTTP_TIMEOUT = 180
MAX_FALLBACK_CYCLES = 24   # 최대 6일 전 run까지 fallback

INTERVAL_RE = re.compile(r":(\d+)-(\d+)\s*hour")


# =========================
# basic helpers
# =========================

def as_dt(x: Any) -> dt:
    if isinstance(x, dt):
        return x
    if isinstance(x, str):
        return dt.strptime(x, "%Y-%m-%d %H:%M:%S")
    raise TypeError(f"Unsupported datetime type: {type(x)} -> {x}")


def get_cfg_datetime(cfg: Dict[str, Any], *keys: str) -> dt:
    for key in keys:
        if key in cfg:
            return as_dt(cfg[key])
    raise KeyError(f"Missing datetime key. Tried: {keys}")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def require_wgrib2() -> None:
    from shutil import which
    if which("wgrib2") is None:
        raise RuntimeError("wgrib2 not found in PATH")


def pick_run_date(cfg: Dict[str, Any], t_st: dt, t_ed: dt) -> str:
    cand_keys = ["initdate", "init_date", "run_date", "RUN_DATE", "RUN_DATE_ISO"]
    for key in cand_keys:
        if key in cfg and isinstance(cfg[key], str):
            s = cfg[key].strip()
            for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%Y%m%d"):
                try:
                    return dt.strptime(s, fmt).strftime("%Y%m%d")
                except Exception:
                    pass
    mid = t_st + (t_ed - t_st) / 2
    return mid.strftime("%Y%m%d")


def make_regular_times(t_st: dt, t_ed: dt, step_hours: int) -> List[dt]:
    step = timedelta(hours=step_hours)
    cur = t_st
    out: List[dt] = []
    while cur <= t_ed:
        out.append(cur)
        cur += step
    return out


def floor_to_run(t: dt) -> dt:
    cands = [t.replace(hour=h, minute=0, second=0, microsecond=0) for h in RUN_HOURS]
    past = [x for x in cands if x <= t]
    if past:
        return max(past)
    return (t - timedelta(days=1)).replace(hour=18, minute=0, second=0, microsecond=0)


def fhr_from(run: dt, valid: dt) -> int:
    return int((valid - run).total_seconds() / 3600)


def gfs_url(run: dt, fhr: int) -> str:
    ymd = run.strftime("%Y%m%d")
    hh = run.strftime("%H")
    return f"{S3_BASE}/gfs.{ymd}/{hh}/atmos/gfs.t{hh}z.pgrb2.{RES}.f{fhr:03d}"


def run_cmd(cmd: List[str], inp: Optional[str] = None) -> str:
    p = subprocess.run(
        cmd,
        input=inp,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout


# =========================
# coords / subset
# =========================

def detect_latlon_names(ds: xr.Dataset) -> Tuple[str, str]:
    for cand in ("lat", "latitude", "Latitude"):
        if cand in ds.coords:
            lat_name = cand
            break
    else:
        raise RuntimeError("latitude coord not found")

    for cand in ("lon", "longitude", "Longitude"):
        if cand in ds.coords:
            lon_name = cand
            break
    else:
        raise RuntimeError("longitude coord not found")

    return lat_name, lon_name


def normalize_lon_range(ds: xr.Dataset, lon_name: str, lon_rng: Tuple[float, float]) -> Tuple[float, float]:
    lon0, lon1 = lon_rng
    lonvals = ds[lon_name].values
    if lonvals.size == 0:
        return lon_rng

    lon_min = float(np.nanmin(lonvals))
    lon_max = float(np.nanmax(lonvals))

    # dataset이 0~360 체계면 음수 lon 입력을 0~360으로 변환
    if lon_min >= 0 and lon_max > 180:
        if lon0 < 0:
            lon0 += 360
        if lon1 < 0:
            lon1 += 360

    return (lon0, lon1)


def safe_subset(ds: xr.Dataset, lat_rng: Tuple[float, float], lon_rng: Tuple[float, float]) -> xr.Dataset:
    lat_name, lon_name = detect_latlon_names(ds)

    if ds[lat_name].size >= 2 and float(ds[lat_name][0]) > float(ds[lat_name][-1]):
        ds = ds.sortby(lat_name)
    if ds[lon_name].size >= 2 and float(ds[lon_name][0]) > float(ds[lon_name][-1]):
        ds = ds.sortby(lon_name)

    lon_rng = normalize_lon_range(ds, lon_name, lon_rng)
    lat0, lat1 = lat_rng
    lon0, lon1 = lon_rng

    return ds.sel({
        lat_name: slice(min(lat0, lat1), max(lat0, lat1)),
        lon_name: slice(min(lon0, lon1), max(lon0, lon1)),
    })


# =========================
# download / fallback
# =========================

def local_raw_path(raw_dir: str, run: dt, fhr: int) -> str:
    return os.path.join(raw_dir, f"{run:%Y%m%d%H}_f{fhr:03d}.grib2")


def try_download(url: str, out_path: str) -> bool:
    if os.path.exists(out_path):
        return True

    tmp_path = out_path + ".part"
    try:
        with requests.get(url, stream=True, timeout=HTTP_TIMEOUT) as r:
            if r.status_code == 404:
                return False
            r.raise_for_status()

            print("download:", url)

            with open(tmp_path, "wb") as f:
                for chunk in r.iter_content(1024 * 1024):
                    if chunk:
                        f.write(chunk)

        os.replace(tmp_path, out_path)
        return True

    finally:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass


def fetch_gfs_file(initial_run: dt, initial_fhr: int, raw_dir: str, label: str) -> Tuple[dt, int, str]:
    """
    planned run/fhr로 먼저 시도하고,
    없으면 run을 6시간씩 뒤로, fhr은 6시간씩 앞으로 보내며 fallback.
    """
    run = initial_run
    fhr = initial_fhr

    for n in range(MAX_FALLBACK_CYCLES + 1):
        if fhr < 0 or fhr > 384:
            break

        raw_path = local_raw_path(raw_dir, run, fhr)

        if os.path.exists(raw_path):
            if n > 0:
                print(f"[INFO] {label} fallback cache -> {run} f{fhr:03d}")
            return run, fhr, raw_path

        url = gfs_url(run, fhr)
        ok = try_download(url, raw_path)

        if ok:
            if n > 0:
                print(f"[INFO] {label} fallback use  -> {run} f{fhr:03d}")
            return run, fhr, raw_path

        print(f"[WARN] {label} missing {run} f{fhr:03d} -> fallback")
        run = run - timedelta(hours=6)
        fhr = fhr + 6

    raise RuntimeError(f"No available GFS file found for {label}: start={initial_run} f{initial_fhr:03d}")


def fetch_exact_file(run: dt, fhr: int, raw_dir: str, label: str) -> str:
    if fhr < 0 or fhr > 384:
        raise RuntimeError(f"{label}: invalid exact fhr={fhr}")
    raw_path = local_raw_path(raw_dir, run, fhr)
    if os.path.exists(raw_path):
        return raw_path

    url = gfs_url(run, fhr)
    ok = try_download(url, raw_path)
    if not ok:
        raise RuntimeError(f"{label}: exact file missing for {run} f{fhr:03d}")
    return raw_path


# =========================
# inventory / selection
# =========================

def inventory_lines(grib_path: str, cache: Dict[str, List[str]]) -> List[str]:
    if grib_path not in cache:
        cache[grib_path] = run_cmd(["wgrib2", grib_path]).splitlines()
    return cache[grib_path]


def select_instant_lines(inv: List[str]) -> List[str]:
    patterns = [
        ":UGRD:10 m above ground:",
        ":VGRD:10 m above ground:",
        ":TMP:2 m above ground:",
        ":SPFH:2 m above ground:",
        ":PRMSL:mean sea level:",
    ]
    out: List[str] = []
    for pat in patterns:
        found = None
        for ln in inv:
            if pat in ln:
                found = ln
                break
        if found is None:
            raise RuntimeError(f"Missing instant variable pattern: {pat}")
        out.append(found)
    return out


def parse_interval(line: str) -> Optional[Tuple[int, int]]:
    m = INTERVAL_RE.search(line)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def find_exact_interval_line(lines: List[str], start: int, end: int) -> Optional[str]:
    for ln in lines:
        iv = parse_interval(ln)
        if iv == (start, end):
            return ln
    return None


def find_covering_interval_line(lines: List[str], target_start: int, target_end: int) -> Optional[Tuple[int, str]]:
    """
    target=(a,b) 를 만들기 위해 end=b 인 더 긴 interval (s,b), s<a 를 찾는다.
    여러 개면 s가 가장 큰 것(가장 가까운 것)을 선택.
    예:
      target 9-12 -> 6-12 선택
      target 3-6  -> 0-6 선택
    """
    cands: List[Tuple[int, str]] = []
    for ln in lines:
        iv = parse_interval(ln)
        if iv is None:
            continue
        s, e = iv
        if e == target_end and s < target_start:
            cands.append((s, ln))

    if not cands:
        return None

    return max(cands, key=lambda x: x[0])


# =========================
# grib -> dataset
# =========================

def extract_lines_to_dataset(
    grib_in: str,
    selected_lines: List[str],
    grib_out: str,
    nc_out: str,
) -> xr.Dataset:
    if not selected_lines:
        raise RuntimeError("No lines selected for extraction")

    if os.path.exists(grib_out):
        os.remove(grib_out)
    if os.path.exists(nc_out):
        os.remove(nc_out)

    run_cmd(["wgrib2", grib_in, "-i", "-grib_out", grib_out], "\n".join(selected_lines) + "\n")
    run_cmd(["wgrib2", grib_out, "-netcdf", nc_out])

    with xr.open_dataset(nc_out) as src:
        ds = src.load()

    return ds


def find_var_by_prefix(ds: xr.Dataset, prefix: str) -> str:
    prefix_u = prefix.upper()
    cands = [v for v in ds.data_vars if v.upper().startswith(prefix_u)]
    if not cands:
        raise RuntimeError(f"Cannot find variable with prefix {prefix} in {list(ds.data_vars)}")
    return cands[0]


def extract_one_var_dataset(
    raw_grib: str,
    line: str,
    raw_dir: str,
    tag: str,
) -> xr.Dataset:
    grib_out = os.path.join(raw_dir, f"{tag}.grib2")
    nc_out = os.path.join(raw_dir, f"{tag}.nc")
    return extract_lines_to_dataset(raw_grib, [line], grib_out, nc_out)


def derive_flux_or_rain(
    ds_start: xr.Dataset,
    ds_end: xr.Dataset,
    var_prefix: str,
    block_start: int,
    mid_time: int,
    end_time: int,
) -> xr.Dataset:
    """
    (block_start, end_time) 와 (block_start, mid_time) 로부터
    (mid_time, end_time) 를 복원한다.

    예:
      0-6 과 0-3 -> 3-6
      6-12 와 6-9 -> 9-12
    """
    vs = find_var_by_prefix(ds_start, var_prefix)
    ve = find_var_by_prefix(ds_end, var_prefix)

    out = xr.Dataset(coords={k: v for k, v in ds_end.coords.items()})

    if var_prefix == "APCP":
        out[ve] = ds_end[ve] - ds_start[vs]
    else:
        dur_end = end_time - block_start
        dur_start = mid_time - block_start
        dur_target = end_time - mid_time
        out[ve] = (ds_end[ve] * dur_end - ds_start[vs] * dur_start) / dur_target

    return out


def build_cumulative_dataset(
    raw_c_end: str,
    act_run_c: dt,
    act_fhr_c: int,
    raw_dir: str,
    inv_cache: Dict[str, List[str]],
) -> xr.Dataset:
    """
    cumulative target = (act_fhr_c-3, act_fhr_c]

    우선 exact a-b 를 찾고,
    없으면 (s-b) 와 (s-a) 조합으로 복원한다.
    """
    a = act_fhr_c - STEP_HOURS
    b = act_fhr_c

    specs = {
        "DSWRF": ":DSWRF:surface:",
        "DLWRF": ":DLWRF:surface:",
        "APCP":  ":APCP:surface:",
    }

    inv_end = inventory_lines(raw_c_end, inv_cache)
    out = xr.Dataset()

    for vprefix, pat in specs.items():
        matches_end = [ln for ln in inv_end if pat in ln]

        # 1) exact a-b
        exact_line = find_exact_interval_line(matches_end, a, b)
        if exact_line is not None:
            ds_var = extract_one_var_dataset(
                raw_c_end,
                exact_line,
                raw_dir,
                f"cum_{vprefix}_{a}_{b}",
            )
            vname = find_var_by_prefix(ds_var, vprefix)
            out[vname] = ds_var[vname]
            print(f"[INFO] cumulative {vprefix}: exact {a}-{b}")
            continue

        # 2) covering (s-b) + start-file의 (s-a) 로 복원
        cover = find_covering_interval_line(matches_end, a, b)
        if cover is None:
            raise RuntimeError(f"Cannot resolve cumulative interval for {vprefix}: need {a}-{b} hour")

        s, end_line = cover

        raw_c_start = fetch_exact_file(act_run_c, a, raw_dir, f"{vprefix}-start")
        inv_start = inventory_lines(raw_c_start, inv_cache)
        matches_start = [ln for ln in inv_start if pat in ln]

        start_line = find_exact_interval_line(matches_start, s, a)
        if start_line is None:
            raise RuntimeError(
                f"Cannot resolve cumulative interval for {vprefix}: "
                f"have {s}-{b} but missing {s}-{a}"
            )

        ds_start = extract_one_var_dataset(
            raw_c_start,
            start_line,
            raw_dir,
            f"cum_{vprefix}_{s}_{a}",
        )
        ds_end = extract_one_var_dataset(
            raw_c_end,
            end_line,
            raw_dir,
            f"cum_{vprefix}_{s}_{b}",
        )

        ds_var = derive_flux_or_rain(ds_start, ds_end, vprefix, s, a, b)
        vname = find_var_by_prefix(ds_var, vprefix)
        out[vname] = ds_var[vname]

        print(f"[INFO] cumulative {vprefix}: derived {a}-{b} from {s}-{a} and {s}-{b}")

    return out


# =========================
# main
# =========================

def main() -> int:
    require_wgrib2()

    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t_st = get_cfg_datetime(cfg, "bry_start_date", "bry_start_time")
    t_ed = get_cfg_datetime(cfg, "bry_end_date", "bry_end_time")

    lat_rng = tuple(cfg["region"]["lat"])
    lon_rng = tuple(cfg["region"]["lon"])

    base_dir = cfg["output"]["base_dir"]
    gfs_dir = os.path.join(base_dir, "gfs")
    raw_dir = os.path.join(gfs_dir, "raw")

    ensure_dir(gfs_dir)
    ensure_dir(raw_dir)

    out_stamp = pick_run_date(cfg, t_st, t_ed)
    out_path = os.path.join(gfs_dir, f"gfs_korea_{out_stamp}.nc")

    print("[INFO] window:", t_st, "→", t_ed)
    print("[INFO] output:", out_path)

    target_times = make_regular_times(t_st, t_ed, STEP_HOURS)
    pieces: List[xr.Dataset] = []
    inv_cache: Dict[str, List[str]] = {}

    for t in target_times:
        plan_run_i = floor_to_run(t)
        plan_fhr_i = fhr_from(plan_run_i, t)

        t_prev = t - timedelta(hours=STEP_HOURS)
        plan_run_c = floor_to_run(t_prev)
        plan_fhr_c = fhr_from(plan_run_c, t)

        print(f"\nTARGET {t}")
        print("instant planned   :", plan_run_i, plan_fhr_i)
        print("cumulative planned:", plan_run_c, plan_fhr_c)

        act_run_i, act_fhr_i, raw_i = fetch_gfs_file(plan_run_i, plan_fhr_i, raw_dir, "instant")
        act_run_c, act_fhr_c, raw_c_end = fetch_gfs_file(plan_run_c, plan_fhr_c, raw_dir, "cumulative")

        print("instant use       :", act_run_i, act_fhr_i)
        print("cumulative use    :", act_run_c, act_fhr_c)

        # ---- instant dataset
        inv_i = inventory_lines(raw_i, inv_cache)
        inst_lines = select_instant_lines(inv_i)
        ds_i = extract_lines_to_dataset(
            raw_i,
            inst_lines,
            os.path.join(raw_dir, "inst.grib2"),
            os.path.join(raw_dir, "inst.nc"),
        )

        # ---- cumulative dataset
        ds_c = build_cumulative_dataset(
            raw_c_end=raw_c_end,
            act_run_c=act_run_c,
            act_fhr_c=act_fhr_c,
            raw_dir=raw_dir,
            inv_cache=inv_cache,
        )

        ds = xr.merge([ds_i, ds_c], compat="override")
        ds = safe_subset(ds, lat_rng, lon_rng)

        days = (t - TIME_REF).total_seconds() / 86400.0

        if ds.sizes.get("time", 0) > 1:
            ds = ds.isel(time=0)

        ds = ds.assign_coords(time=("time", [days]))
        ds["time"].attrs["units"] = "days since 2000-01-01 00:00:00"
        ds["time"].attrs["calendar"] = "proleptic_gregorian"

        pieces.append(ds)

    final = xr.concat(pieces, dim="time").sortby("time")

    tmp_out = out_path + ".tmp"
    final.to_netcdf(tmp_out)
    os.replace(tmp_out, out_path)

    print("\n[OK] saved:", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
