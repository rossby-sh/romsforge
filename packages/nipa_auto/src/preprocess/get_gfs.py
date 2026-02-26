
#!/usr/bin/env python3
from __future__ import annotations

import os
import time
import yaml
import shutil
import subprocess
from datetime import datetime as dt, timedelta
from typing import Any, List, Tuple, Optional

import requests


# =========================
# utils
# =========================
def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def as_dt(x: Any) -> dt:
    if isinstance(x, dt):
        return x
    if isinstance(x, str):
        return dt.strptime(x, "%Y-%m-%d %H:%M:%S")
    raise TypeError(f"Unsupported time type: {type(x)} -> {x}")

def retry(fn, tries: int = 5, wait_sec: int = 120, name: str = "task"):
    last_err = None
    for n in range(1, tries + 1):
        try:
            return fn()
        except Exception as e:
            last_err = e
            if n == tries:
                raise
            print(f"[WARN] {name} failed ({n}/{tries}): {e} -> sleep {wait_sec}s", flush=True)
            time.sleep(wait_sec)
    raise last_err  # unreachable

def download_file(url: str, out_path: str, timeout: int = 180) -> None:
    tmp = out_path + ".part"
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    os.replace(tmp, out_path)

def run(cmd: List[str], *, inp: Optional[str] = None) -> str:
    p = subprocess.run(
        cmd,
        input=inp,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout

def require_wgrib2():
    if shutil.which("wgrib2") is None:
        raise RuntimeError("wgrib2 not found in PATH.")


# =========================
# GFS (NOMADS) run selection
# =========================
GFS_RUN_HOURS = (0, 6, 12, 18)

def pick_latest_run_at_or_before(cur: dt) -> dt:
    candidates = [cur.replace(hour=hh, minute=0, second=0, microsecond=0) for hh in GFS_RUN_HOURS]
    past = [rt for rt in candidates if rt <= cur]
    if past:
        return max(past)
    return (cur - timedelta(days=1)).replace(hour=18, minute=0, second=0, microsecond=0)

def build_gfs_url_for_run(run_time: dt, fhr: int) -> str:
    run_date = run_time.strftime("%Y%m%d")
    run_hh = run_time.strftime("%H")
    return (
        "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/"
        f"gfs.{run_date}/{run_hh}/atmos/"
        f"gfs.t{run_hh}z.pgrb2.0p50.f{fhr:03d}"
    )


# =========================
# wgrib2 helpers
# =========================
def wgrib2_lines(grib_path: str) -> List[str]:
    out = run(["wgrib2", grib_path])
    return out.splitlines()

def grep_lines(lines: List[str], pattern: str) -> List[str]:
    return [ln for ln in lines if pattern in ln]

def pick_first(lines: List[str]) -> Optional[str]:
    return lines[0] if lines else None

def pick_interval_line(lines: List[str], fhr: int) -> Optional[str]:
    """
    pick line that matches (fhr-3)-(fhr) hour interval.
    For cumulative vars we use fhr=3 => "0-3 hour".
    """
    if not lines:
        return None
    if fhr <= 0:
        return lines[0]
    a = fhr - 3
    b = fhr
    key = f":{a}-{b} hour"
    for ln in lines:
        if key in ln:
            return ln
    return lines[0]

def wgrib2_extract_records(grib_path: str, out_grib: str, match_lines: List[str]) -> None:
    if not match_lines:
        raise RuntimeError("No match lines to extract.")
    inp = "\n".join(match_lines) + "\n"
    run(["wgrib2", grib_path, "-i", "-grib_out", out_grib], inp=inp)

def wgrib2_append_records(grib_path: str, out_grib: str, match_lines: List[str]) -> None:
    if not match_lines:
        return
    inp = "\n".join(match_lines) + "\n"
    run(["wgrib2", grib_path, "-i", "-append", "-grib_out", out_grib], inp=inp)

def wgrib2_to_netcdf(grib_path: str, nc_path: str) -> None:
    run(["wgrib2", grib_path, "-netcdf", nc_path])


# =========================
# selection rules
# =========================
def select_instant_lines(inv_f000: List[str]) -> Tuple[dict, List[str]]:
    """
    instant variables from f000 file
    """
    pat_u10   = ":UGRD:10 m above ground:"
    pat_v10   = ":VGRD:10 m above ground:"
    pat_tair  = ":TMP:2 m above ground:"
    pat_qair  = ":SPFH:2 m above ground:"
    pat_pair  = ":PRMSL:mean sea level:"

    # cloud: prefer total cloud cover entire atmosphere
    pat_cloud_best = ":TCDC:entire atmosphere:"
    pat_cloud_any  = ":TCDC:"

    def one(pat: str) -> Optional[str]:
        lines = grep_lines(inv_f000, pat)
        return lines[0] if lines else None

    sel = {}
    sel["u10"]  = one(pat_u10)
    sel["v10"]  = one(pat_v10)
    sel["tair"] = one(pat_tair)
    sel["qair"] = one(pat_qair)
    sel["pair"] = one(pat_pair)

    cloud = one(pat_cloud_best)
    if cloud is None:
        cloud = one(pat_cloud_any)
    sel["cloud"] = cloud  # optional

    required = ["u10", "v10", "tair", "qair", "pair"]
    missing = [k for k in required if sel.get(k) is None]
    return sel, missing

def select_cumulative_lines(inv_f003: List[str]) -> Tuple[dict, List[str]]:
    """
    cumulative/interval variables from f003 file (0-3 hour interval)
    """
    pat_swrad = ":DSWRF:surface:"
    pat_lwrad = ":DLWRF:surface:"
    pat_rain  = ":APCP:surface:"

    def interval(pat: str) -> Optional[str]:
        lines = grep_lines(inv_f003, pat)
        return pick_interval_line(lines, fhr=3) if lines else None

    sel = {}
    sel["swrad"] = interval(pat_swrad)
    sel["lwrad"] = interval(pat_lwrad)
    sel["rain"]  = interval(pat_rain)

    required = ["swrad", "lwrad", "rain"]
    missing = [k for k in required if sel.get(k) is None]
    return sel, missing


# =========================
# postprocess nc: collapse time=2 -> time=1
# =========================
def collapse_to_single_time(nc_path: str, run_time: dt) -> None:
    """
    wgrib2 netcdf의 time=2 (f000/f003) 문제를 해결:
      - instant(UGRD/VGRD/TMP/SPFH/PRMSL/TCDC)은 time=0
      - cumulative(DSWRF/DLWRF/APCP)은 time=1
    로 슬라이스해서 time=1개짜리로 다시 저장한다.
    """
    import numpy as np
    import xarray as xr

    ds = xr.open_dataset(nc_path)

    if "time" not in ds.dims:
        ds.close()
        return
    tlen = ds.sizes.get("time", 0)
    if tlen <= 1:
        ds.close()
        return

    inst_keys = ("UGRD", "VGRD", "TMP", "SPFH", "PRMSL", "TCDC")
    cum_keys  = ("DSWRF", "DLWRF", "APCP")

    inst_vars, cum_vars, other_vars = [], [], []
    for v in ds.data_vars:
        name = v.upper()
        if any(k in name for k in cum_keys):
            cum_vars.append(v)
        elif any(k in name for k in inst_keys):
            inst_vars.append(v)
        else:
            other_vars.append(v)

    # time index: assume [0]=f000, [1]=f003 after append
    ds_inst = ds[inst_vars].isel(time=0) if inst_vars else None
    ds_cum  = ds[cum_vars].isel(time=1)  if cum_vars  else None
    ds_oth  = ds[other_vars].isel(time=0) if other_vars else None

    parts = []
    for part in (ds_inst, ds_cum, ds_oth):
        if part is None:
            continue
        part = part.expand_dims(time=[np.datetime64(run_time)])
        parts.append(part)

    if not parts:
        ds.close()
        return

    out = xr.merge(parts, compat="override")
    out["time"].attrs["long_name"] = "time"

    tmp = nc_path + ".tmp"
    out.to_netcdf(tmp)
    out.close()
    ds.close()
    os.replace(tmp, nc_path)


# =========================
# main
# =========================
def main() -> int:
    require_wgrib2()

    HERE = os.path.dirname(os.path.abspath(__file__))
    CFG_PATH = os.path.join(HERE, "config.yaml")

    with open(CFG_PATH) as f:
        cfg = yaml.safe_load(f)

    t_st = as_dt(cfg["bry_start_date"])
    t_ed = as_dt(cfg["bry_end_date"])

    base_dir = cfg["output"]["base_dir"]
    out_dir = os.path.join(base_dir, "gfs")
    raw_dir = os.path.join(out_dir, "raw")
    thin_dir = os.path.join(out_dir, "thin")
    ensure_dir(out_dir)
    ensure_dir(raw_dir)
    ensure_dir(thin_dir)

    print("[INFO] config:", CFG_PATH)
    print("[INFO] bry window:", t_st, "to", t_ed)
    print("[INFO] output:", out_dir)
    print("[INFO] step: 6 hours")
    print("[INFO] policy: instant from f000, cumulative from f003 (same run)")
    print("[INFO] note: cumulative time is collapsed to f000 run_time in final nc")

    step = timedelta(hours=6)
    cur = t_st

    while cur <= t_ed:
        run_time = pick_latest_run_at_or_before(cur)
        run_stamp = run_time.strftime("%Y%m%d%H")

        grib_f000 = os.path.join(raw_dir,  f"gfs_run_{run_stamp}_f000.grib2")
        grib_f003 = os.path.join(raw_dir,  f"gfs_run_{run_stamp}_f003.grib2")

        grib_thin = os.path.join(thin_dir, f"gfs_run_{run_stamp}_thin.grib2")
        nc_out    = os.path.join(out_dir,  f"gfs_run_{run_stamp}.nc")

        # ---- download f000 ----
        if not os.path.exists(grib_f000):
            url0 = build_gfs_url_for_run(run_time, fhr=0)
            print(f"\n[INFO] run={run_stamp} download f000")
            print("[INFO] url:", url0)

            def _dl0():
                download_file(url0, grib_f000, timeout=180)

            try:
                retry(_dl0, tries=5, wait_sec=120, name=f"download {run_stamp} f000")
                if os.path.getsize(grib_f000) < 1024 * 200:
                    raise RuntimeError("Downloaded f000 too small (likely error response).")
                print("[OK] raw saved:", grib_f000)
            except Exception as e:
                print("[ERR] download f000 failed:", run_stamp, e)
                cur += step
                continue
        else:
            print(f"\n[SKIP] raw f000 exists: {grib_f000}")

        # ---- download f003 ----
        if not os.path.exists(grib_f003):
            url3 = build_gfs_url_for_run(run_time, fhr=3)
            print(f"[INFO] run={run_stamp} download f003")
            print("[INFO] url:", url3)

            def _dl3():
                download_file(url3, grib_f003, timeout=180)

            try:
                retry(_dl3, tries=5, wait_sec=120, name=f"download {run_stamp} f003")
                if os.path.getsize(grib_f003) < 1024 * 200:
                    raise RuntimeError("Downloaded f003 too small (likely error response).")
                print("[OK] raw saved:", grib_f003)
            except Exception as e:
                print("[ERR] download f003 failed:", run_stamp, e)
                cur += step
                continue
        else:
            print(f"[SKIP] raw f003 exists: {grib_f003}")

        # ---- build thin (f000 then append f003) ----
        if not os.path.exists(grib_thin):
            inv0 = wgrib2_lines(grib_f000)
            inv3 = wgrib2_lines(grib_f003)

            sel0, miss0 = select_instant_lines(inv0)
            sel3, miss3 = select_cumulative_lines(inv3)

            if miss0:
                print("[ERR] missing required instant from f000:", miss0)
                print("[ERR] skip run:", run_stamp)
                cur += step
                continue
            if miss3:
                print("[ERR] missing required cumulative from f003:", miss3)
                print("[ERR] skip run:", run_stamp)
                cur += step
                continue

            # f000 -> thin
            lines0 = [sel0[k] for k in ["u10", "v10", "tair", "qair", "pair"] if sel0.get(k)]
            if sel0.get("cloud"):
                lines0.append(sel0["cloud"])

            try:
                wgrib2_extract_records(grib_f000, grib_thin, lines0)
            except Exception as e:
                print("[ERR] thin build (f000 part) failed:", run_stamp, e)
                cur += step
                continue

            # f003 append
            lines3 = [sel3[k] for k in ["swrad", "lwrad", "rain"] if sel3.get(k)]
            try:
                wgrib2_append_records(grib_f003, grib_thin, lines3)
            except Exception as e:
                print("[ERR] thin build (f003 part) failed:", run_stamp, e)
                cur += step
                continue

            print("[OK] thin grib saved:", grib_thin)
        else:
            print(f"[SKIP] thin exists: {grib_thin}")

        # ---- thin -> netcdf ----
        if not os.path.exists(nc_out):
            try:
                wgrib2_to_netcdf(grib_thin, nc_out)
                print("[OK] nc saved:", nc_out)
            except Exception as e:
                print("[ERR] netcdf convert failed:", run_stamp, e)
                cur += step
                continue
        else:
            print(f"[SKIP] nc exists: {nc_out}")

        # ---- collapse time dimension to run_time (force single time) ----
        try:
            collapse_to_single_time(nc_out, run_time)
            print("[OK] nc time collapsed to single run_time:", nc_out)
        except Exception as e:
            print("[WARN] nc time-collapse failed (kept original):", e)

        cur += step

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
