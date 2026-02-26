
#!/usr/bin/env python3
# get_gfs2.py
#
# 목표:
# - window: 2026-01-22 21:00:00 ~ 2026-01-24 03:00:00 (3시간 간격)
# - 각 valid time마다 최신 GFS cycle(00/06/12/18) 중 fhr<=30으로 커버 가능한 run 선택
# - 변수:
#   instant: u10,v10,t2,q2,prmsl,cloud
#   interval(3h): swrad(DSWRF), lwrad(DLWRF), rain(APCP)
# - fhr=0일 때 누적/구간 변수가 없으면:
#   prev run = run_time - 6h, prev_fhr=6 (3-6 hour interval = vt-3h ~ vt)에서 가져옴
#
# 필요:
# - wgrib2
# - python: requests, xarray, numpy, netCDF4 (xarray backend)

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime as dt, timedelta
from typing import Optional, List, Dict, Tuple

import requests


# =========================
# config
# =========================
START = "2026-01-22 21:00:00"
END   = "2026-01-24 03:00:00"
STEP_HOURS = 3

MAX_FHR = 30
GRID = "0p50"

OUT_DIR = "./gfs_out"
OUT_NC  = os.path.join(OUT_DIR, "gfs_20260122_2100__20260124_0300.nc")

NOMADS_BASE = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
GFS_RUN_HOURS = (0, 6, 12, 18)


# =========================
# utils
# =========================
def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def as_dt(s: str) -> dt:
    return dt.strptime(s, "%Y-%m-%d %H:%M:%S")

def require_cmd(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"{name} not found in PATH.")

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

def download_file(url: str, out_path: str, timeout: int = 180) -> None:
    tmp = out_path + ".part"
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    os.replace(tmp, out_path)

def is_grib2(path: str) -> bool:
    try:
        with open(path, "rb") as f:
            head = f.read(4)
        return head == b"GRIB"
    except Exception:
        return False


# =========================
# GFS pick / URL
# =========================
def build_gfs_url(run_time: dt, fhr: int) -> str:
    run_date = run_time.strftime("%Y%m%d")
    run_hh = run_time.strftime("%H")
    return (
        f"{NOMADS_BASE}/gfs.{run_date}/{run_hh}/atmos/"
        f"gfs.t{run_hh}z.pgrb2.{GRID}.f{fhr:03d}"
    )

@dataclass(frozen=True)
class RunPick:
    run_time: dt
    fhr: int

def pick_run_for_valid_time(valid_time: dt, max_fhr: int) -> RunPick:
    """
    valid_time을 커버하는 최신 run_time(00/06/12/18)과 fhr 선택.
    조건:
      - run_time <= valid_time
      - fhr = (valid_time - run_time) hours
      - fhr % 3 == 0
      - 0 <= fhr <= max_fhr
    """
    candidates: List[RunPick] = []

    # max_fhr=30이면 하루/전날 정도면 충분
    for day_offset in (0, 1):
        day = (valid_time - timedelta(days=day_offset)).date()
        for hh in GFS_RUN_HOURS:
            rt = dt(day.year, day.month, day.day, hh, 0, 0)
            if rt > valid_time:
                continue
            fhr = int((valid_time - rt).total_seconds() // 3600)
            if fhr < 0 or fhr > max_fhr:
                continue
            if fhr % 3 != 0:
                continue
            candidates.append(RunPick(rt, fhr))

    if not candidates:
        raise RuntimeError(
            f"cannot pick run for valid_time={valid_time} within fhr<= {max_fhr} and 3h step"
        )

    # 최신 run_time 우선
    best = max(candidates, key=lambda x: x.run_time)
    return best


# =========================
# wgrib2 helpers
# =========================
def wgrib2_lines(grib_path: str) -> List[str]:
    out = run(["wgrib2", grib_path])
    return out.splitlines()

def grep_lines(lines: List[str], pattern: str) -> List[str]:
    return [ln for ln in lines if pattern in ln]

def pick_interval_line(lines: List[str], fhr: int) -> Optional[str]:
    """
    interval/cumulative vars: pick "{fhr-3}-{fhr} hour" line.
    """
    if not lines:
        return None
    if fhr <= 0:
        return None
    a = fhr - 3
    b = fhr
    key = f":{a}-{b} hour"
    for ln in lines:
        if key in ln:
            return ln
    return None

def wgrib2_extract(grib_path: str, out_grib: str, match_lines: List[str], *, append: bool = False) -> None:
    """
    match_lines (wgrib2 inventory line strings)만 골라 grib2로 저장.
    append=True면 out_grib 뒤에 append.
    """
    if not match_lines:
        raise RuntimeError("No match lines to extract.")
    inp = "\n".join(match_lines) + "\n"
    cmd = ["wgrib2", grib_path, "-i"]
    if append:
        cmd += ["-append"]
    cmd += ["-grib_out", out_grib]
    run(cmd, inp=inp)

def wgrib2_to_netcdf(grib_path: str, nc_path: str) -> None:
    run(["wgrib2", grib_path, "-netcdf", nc_path])


# =========================
# selection rules
# =========================
def select_instant(inv: List[str]) -> Tuple[Dict[str, str], List[str]]:
    pat_u10  = ":UGRD:10 m above ground:"
    pat_v10  = ":VGRD:10 m above ground:"
    pat_tair = ":TMP:2 m above ground:"
    pat_qair = ":SPFH:2 m above ground:"
    pat_pair = ":PRMSL:mean sea level:"

    pat_cloud_best = ":TCDC:entire atmosphere:"
    pat_cloud_any  = ":TCDC:"

    def first(pat: str) -> Optional[str]:
        m = grep_lines(inv, pat)
        return m[0] if m else None

    sel: Dict[str, Optional[str]] = {}
    sel["u10"]  = first(pat_u10)
    sel["v10"]  = first(pat_v10)
    sel["tair"] = first(pat_tair)
    sel["qair"] = first(pat_qair)
    sel["pair"] = first(pat_pair)

    cloud = first(pat_cloud_best)
    if cloud is None:
        cloud = first(pat_cloud_any)
    sel["cloud"] = cloud

    required = ["u10", "v10", "tair", "qair", "pair"]
    missing = [k for k in required if sel.get(k) is None]

    out: Dict[str, str] = {k: v for k, v in sel.items() if v is not None}
    return out, missing

def select_cumulative(inv: List[str], fhr: int) -> Tuple[Dict[str, str], List[str]]:
    pat_swrad = ":DSWRF:surface:"
    pat_lwrad = ":DLWRF:surface:"
    pat_rain  = ":APCP:surface:"

    def interval(pat: str) -> Optional[str]:
        m = grep_lines(inv, pat)
        return pick_interval_line(m, fhr=fhr)

    if fhr <= 0:
        return {}, ["swrad", "lwrad", "rain"]

    sel: Dict[str, Optional[str]] = {}
    sel["swrad"] = interval(pat_swrad)
    sel["lwrad"] = interval(pat_lwrad)
    sel["rain"]  = interval(pat_rain)

    required = ["swrad", "lwrad", "rain"]
    missing = [k for k in required if sel.get(k) is None]
    out: Dict[str, str] = {k: v for k, v in sel.items() if v is not None}
    return out, missing


# =========================
# netcdf standardization
# =========================
def standardize_step_nc(nc_path: str, valid_time: dt):
    """
    wgrib2 -netcdf 결과에서 변수를 표준 이름으로 rename하고,
    time 좌표를 valid_time으로 지정한다.
    """
    import numpy as np
    import xarray as xr

    ds = xr.open_dataset(nc_path)

    def find_var(keys: Tuple[str, ...]) -> Optional[str]:
        for v in ds.data_vars:
            name = v.upper()
            if all(k in name for k in keys):
                return v
        return None

    mapping = {
        "u10":   find_var(("UGRD", "10")),
        "v10":   find_var(("VGRD", "10")),
        "tair":  find_var(("TMP", "2")),
        "qair":  find_var(("SPFH", "2")),
        "pair":  find_var(("PRMSL",)),
        "cloud": find_var(("TCDC",)),
        "swrad": find_var(("DSWRF",)),
        "lwrad": find_var(("DLWRF",)),
        "rain":  find_var(("APCP",)),
    }

    required = ("u10", "v10", "tair", "qair", "pair", "swrad", "lwrad", "rain")
    missing = [k for k in required if mapping.get(k) is None]
    if missing:
        ds.close()
        raise RuntimeError(f"missing vars in {nc_path}: {missing}")

    keep_vars = [v for v in mapping.values() if v is not None]
    out = ds[keep_vars].copy()

    rename = {mapping[k]: k for k in mapping if mapping[k] is not None}
    out = out.rename(rename)

    # wgrib2 netcdf에는 이미 time dim이 있는 경우가 흔함 -> expand_dims 금지
    out = out.assign_coords(time=[np.datetime64(valid_time)])
    out["time"].attrs["long_name"] = "time"

    ds.close()
    return out


# =========================
# main
# =========================
def main() -> int:
    require_cmd("wgrib2")
    ensure_dir(OUT_DIR)

    t_st = as_dt(START)
    t_ed = as_dt(END)

    raw_dir  = os.path.join(OUT_DIR, "raw")
    step_dir = os.path.join(OUT_DIR, "step")
    ensure_dir(raw_dir)
    ensure_dir(step_dir)

    times: List[dt] = []
    cur = t_st
    while cur <= t_ed:
        times.append(cur)
        cur += timedelta(hours=STEP_HOURS)

    print(f"[INFO] window: {t_st} ~ {t_ed} (step={STEP_HOURS}h) -> n={len(times)}")
    print(f"[INFO] policy: pick latest run covering each valid time within fhr<= {MAX_FHR}")
    print(f"[INFO] out: {OUT_NC}")

    step_datasets = []

    for vt in times:
        pick = pick_run_for_valid_time(vt, max_fhr=MAX_FHR)
        run_stamp = pick.run_time.strftime("%Y%m%d%H")
        fhr = pick.fhr

        grib_path = os.path.join(raw_dir, f"gfs_{run_stamp}_f{fhr:03d}.grib2")
        thin_grib = os.path.join(step_dir, f"thin_{vt.strftime('%Y%m%d%H')}.grib2")
        step_nc   = os.path.join(step_dir, f"step_{vt.strftime('%Y%m%d%H')}.nc")

        url = build_gfs_url(pick.run_time, fhr=fhr)
        print(f"\n[INFO] vt={vt} -> run={run_stamp} fhr={fhr:03d}")
        print(f"[INFO] url={url}")

        # ---- download base grib ----
        if not os.path.exists(grib_path):
            download_file(url, grib_path, timeout=180)
            if (os.path.getsize(grib_path) < 1024 * 200) or (not is_grib2(grib_path)):
                raise RuntimeError(f"download looks invalid: {grib_path}")
            print(f"[OK] raw: {grib_path}")
        else:
            print(f"[SKIP] raw exists: {grib_path}")

        inv = wgrib2_lines(grib_path)

        # ---- instant selection ----
        inst_sel, inst_missing = select_instant(inv)
        if inst_missing:
            raise RuntimeError(f"missing instant in {grib_path}: {inst_missing}")

        use_prev_cum = (fhr == 0)

        # ---- cumulative selection ----
        if not use_prev_cum:
            cum_sel, cum_missing = select_cumulative(inv, fhr=fhr)
            if cum_missing:
                raise RuntimeError(
                    f"missing cumulative for vt={vt} (run={run_stamp} fhr={fhr}): {cum_missing}\n"
                    f"hint: cumulative needs fhr>=3 and '{fhr-3}-{fhr} hour' records"
                )
            prev_grib_path = None
        else:
            # prev run = run_time - 6h, use f006 (3-6 hour interval -> vt-3h ~ vt)
            prev_run_time = pick.run_time - timedelta(hours=6)
            prev_fhr = 6

            prev_stamp = prev_run_time.strftime("%Y%m%d%H")
            prev_grib_path = os.path.join(raw_dir, f"gfs_{prev_stamp}_f{prev_fhr:03d}.grib2")
            prev_url = build_gfs_url(prev_run_time, fhr=prev_fhr)

            print(f"[INFO] cumulative fallback: vt={vt} -> prev run={prev_stamp} fhr={prev_fhr:03d}")
            print(f"[INFO] prev url={prev_url}")

            if not os.path.exists(prev_grib_path):
                download_file(prev_url, prev_grib_path, timeout=180)
                if (os.path.getsize(prev_grib_path) < 1024 * 200) or (not is_grib2(prev_grib_path)):
                    raise RuntimeError(f"download looks invalid: {prev_grib_path}")
                print(f"[OK] prev raw: {prev_grib_path}")
            else:
                print(f"[SKIP] prev raw exists: {prev_grib_path}")

            prev_inv = wgrib2_lines(prev_grib_path)
            cum_sel, cum_missing = select_cumulative(prev_inv, fhr=prev_fhr)
            if cum_missing:
                raise RuntimeError(
                    f"missing cumulative even in prev run for vt={vt} (prev run={prev_stamp} fhr={prev_fhr}): {cum_missing}"
                )

        # ---- build thin grib (inst + cum append) ----
        if not os.path.exists(thin_grib):
            inst_lines = []
            for k in ("u10", "v10", "tair", "qair", "pair"):
                inst_lines.append(inst_sel[k])
            if "cloud" in inst_sel:
                inst_lines.append(inst_sel["cloud"])

            # write inst first
            wgrib2_extract(grib_path, thin_grib, inst_lines, append=False)

            # append cumulative
            cum_lines = [cum_sel[k] for k in ("swrad", "lwrad", "rain")]

            if use_prev_cum:
                assert prev_grib_path is not None
                wgrib2_extract(prev_grib_path, thin_grib, cum_lines, append=True)
            else:
                wgrib2_extract(grib_path, thin_grib, cum_lines, append=True)

            print(f"[OK] thin: {thin_grib}")
        else:
            print(f"[SKIP] thin exists: {thin_grib}")

        # ---- thin -> netcdf ----
        if not os.path.exists(step_nc):
            wgrib2_to_netcdf(thin_grib, step_nc)
            print(f"[OK] step nc: {step_nc}")
        else:
            print(f"[SKIP] step nc exists: {step_nc}")

        # ---- standardize & collect ----
        ds_step = standardize_step_nc(step_nc, valid_time=vt)
        step_datasets.append(ds_step)

    # ---- concat all times ----
    import xarray as xr

    ds_all = xr.concat(step_datasets, dim="time")
    ds_all = ds_all.sortby("time")

    tmp = OUT_NC + ".tmp"
    ds_all.to_netcdf(tmp)
    ds_all.close()

    for d in step_datasets:
        d.close()

    os.replace(tmp, OUT_NC)
    print(f"\n[OK] final nc: {OUT_NC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
