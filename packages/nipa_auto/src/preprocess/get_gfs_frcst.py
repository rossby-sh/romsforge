
#!/usr/bin/env python3
from __future__ import annotations

import os
import time
import yaml
import shutil
import subprocess
from datetime import datetime as dt, timedelta
from typing import Any, List, Tuple, Optional, Dict

import requests
import numpy as np
import xarray as xr


# =========================
# CONFIG (policy)
# =========================
S3_BASE = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
# 해상도는 사이트/시점마다 다를 수 있어서 둘 다 시도
RES_CANDIDATES = ["0p50", "0p25"]

# target time grid
DT_HOURS = 3
WINDOW_START_OFFSET_H = -6
WINDOW_END_OFFSET_H   = 30

# retry
DL_TRIES = 6
DL_WAIT_SEC = 300  # 5 min


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

def require_wgrib2():
    if shutil.which("wgrib2") is None:
        raise RuntimeError("wgrib2 not found in PATH.")

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

def retry(fn, tries: int, wait_sec: int, name: str):
    last = None
    for i in range(1, tries + 1):
        try:
            return fn()
        except Exception as e:
            last = e
            if i == tries:
                raise
            print(f"[WARN] {name} failed ({i}/{tries}): {e} -> sleep {wait_sec}s", flush=True)
            time.sleep(wait_sec)
    raise last  # unreachable

def pick_run_date(cfg: Dict[str, Any], center_dt: dt) -> str:
    # 네 config가 RUN_DATE 같은 걸 들고 있으면 그걸 우선
    for k in ["RUN_DATE", "RUN_DATE_ISO", "run_date", "initdate", "init_date"]:
        if k in cfg and isinstance(cfg[k], str):
            s = cfg[k].strip()
            for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%Y%m%d"):
                try:
                    return dt.strptime(s, fmt).strftime("%Y%m%d")
                except Exception:
                    pass
    return center_dt.strftime("%Y%m%d")


# =========================
# run/time logic
# =========================
RUN_HOURS = (0, 6, 12, 18)

def floor_to_run(t: dt) -> dt:
    """가장 최근 run(00/06/12/18)로 내림"""
    cands = [t.replace(hour=h, minute=0, second=0, microsecond=0) for h in RUN_HOURS]
    past = [x for x in cands if x <= t]
    if past:
        return max(past)
    # 전날 18z
    return (t - timedelta(days=1)).replace(hour=18, minute=0, second=0, microsecond=0)

def fhr_from(run_time: dt, valid: dt) -> int:
    dh = int((valid - run_time).total_seconds() // 3600)
    return dh

def build_s3_url(run_time: dt, fhr: int, res: str) -> str:
    ymd = run_time.strftime("%Y%m%d")
    hh  = run_time.strftime("%H")
    return f"{S3_BASE}/gfs.{ymd}/{hh}/atmos/gfs.t{hh}z.pgrb2.{res}.f{fhr:03d}"

def candidate_urls(run_time: dt, fhr: int) -> List[str]:
    return [build_s3_url(run_time, fhr, res) for res in RES_CANDIDATES]

def download_any(urls: List[str], out_path: str, timeout: int = 180) -> str:
    """
    여러 URL 중 성공하는 첫번째를 다운.
    return: 성공한 url
    """
    headers = {
        "User-Agent": "Mozilla/5.0 (compatible; nipa-auto/1.0)",
        "Accept": "*/*",
        "Connection": "close",
    }
    tmp = out_path + ".part"
    last_err: Optional[Exception] = None

    for url in urls:
        try:
            with requests.get(url, stream=True, timeout=timeout, headers=headers) as r:
                if r.status_code in (403, 404):
                    raise RuntimeError(f"{r.status_code} {r.reason} for url: {url}")
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
            os.replace(tmp, out_path)
            return url
        except Exception as e:
            last_err = e
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass
            continue

    raise last_err if last_err else RuntimeError("download failed")


# =========================
# wgrib2 pick + thin
# =========================
def wgrib2_lines(grib_path: str) -> List[str]:
    return run(["wgrib2", grib_path]).splitlines()

def grep_lines(lines: List[str], pattern: str) -> List[str]:
    return [ln for ln in lines if pattern in ln]

def pick_interval(lines: List[str], fhr: int) -> Optional[str]:
    """
    DSWRF/DLWRF/APCP 같은 interval/accum은 보통
      0-3, 3-6, 6-9 ...
    따라서 fhr=6이면 3-6을 골라야 함.
    """
    if not lines:
        return None
    if fhr < 3:
        return None
    a = fhr - 3
    b = fhr
    key = f":{a}-{b} hour"
    for ln in lines:
        if key in ln:
            return ln
    # fallback: 일단 첫번째
    return lines[0]

def wgrib2_extract(grib_in: str, grib_out: str, selected_lines: List[str]) -> None:
    if not selected_lines:
        raise RuntimeError("No lines selected for extraction.")
    inp = "\n".join(selected_lines) + "\n"
    run(["wgrib2", grib_in, "-i", "-grib_out", grib_out], inp=inp)

def wgrib2_append(grib_in: str, grib_out: str, selected_lines: List[str]) -> None:
    if not selected_lines:
        return
    inp = "\n".join(selected_lines) + "\n"
    run(["wgrib2", grib_in, "-i", "-append", "-grib_out", grib_out], inp=inp)

def wgrib2_to_netcdf(grib_in: str, nc_out: str) -> None:
    run(["wgrib2", grib_in, "-netcdf", nc_out])


def select_instant(inv: List[str]) -> Tuple[List[str], List[str]]:
    """
    instant: u10 v10 t2m spfh2m prmsl cloud
    """
    want = {
        "u10":   ":UGRD:10 m above ground:",
        "v10":   ":VGRD:10 m above ground:",
        "tair":  ":TMP:2 m above ground:",
        "qair":  ":SPFH:2 m above ground:",
        "pair":  ":PRMSL:mean sea level:",
        "cloud": ":TCDC:entire atmosphere:",
    }

    picked: Dict[str, Optional[str]] = {k: None for k in want}
    for k, pat in want.items():
        lst = grep_lines(inv, pat)
        if lst:
            picked[k] = lst[0]

    # cloud는 optional
    missing = [k for k in ["u10","v10","tair","qair","pair"] if picked[k] is None]
    out_lines = [picked[k] for k in ["u10","v10","tair","qair","pair"] if picked[k] is not None]
    if picked["cloud"] is not None:
        out_lines.append(picked["cloud"])
    return out_lines, missing

def select_cumulative(inv: List[str], fhr: int) -> Tuple[List[str], List[str]]:
    """
    cumulative/interval: swrad lwrad rain
    DSWRF/DLWRF: ave, APCP: acc
    """
    pats = {
        "swrad": ":DSWRF:surface:",
        "lwrad": ":DLWRF:surface:",
        "rain":  ":APCP:surface:",
    }

    picked: Dict[str, Optional[str]] = {k: None for k in pats}
    for k, pat in pats.items():
        lst = grep_lines(inv, pat)
        picked[k] = pick_interval(lst, fhr)

    missing = [k for k in ["swrad","lwrad","rain"] if picked[k] is None]
    out_lines = [picked[k] for k in ["swrad","lwrad","rain"] if picked[k] is not None]
    return out_lines, missing


# =========================
# per-time netcdf post
# =========================

def force_single_time_and_set(ds: xr.Dataset, valid_time: dt) -> xr.Dataset:
    """
    최종적으로 time=1개로 만들고, time 좌표를 valid_time 하나로 강제.
    - time이 이미 있으면: isel(time=0)로 1개로 만든 뒤 coords 덮어씀
    - time이 없으면: expand_dims로 생성
    """
    # drop coords we don't want
    for c in ("step", "valid_time"):
        if c in ds.coords:
            ds = ds.drop_vars(c)

    vt = np.datetime64(valid_time)

    if "time" in ds.dims:
        # time을 1개로 축소
        if ds.sizes.get("time", 0) > 1:
            ds = ds.isel(time=0)
        # time 좌표값을 valid_time으로 덮어씀 (차원 유지)
        ds = ds.assign_coords(time=("time", [vt]))
    else:
        ds = ds.expand_dims(time=[vt])

    return ds

# =========================
# main
# =========================
def main() -> int:
    require_wgrib2()

    HERE = os.path.dirname(os.path.abspath(__file__))
    CFG_PATH = os.path.join(HERE, "config.yaml")

    with open(CFG_PATH) as f:
        cfg = yaml.safe_load(f)

    # config에는 이미 UTC 기반 bry_start/end가 들어있다고 가정
    # (너가 main.sh에서 TZ=UTC로 맞추고 만들었으니까)
    bry_st = as_dt(cfg["bry_start_date"])
    bry_ed = as_dt(cfg["bry_end_date"])

    # “00시 기준” center를 bry window에서 추정하지 말고,
    # run_date를 기준으로 00:00을 잡음:
    #   - cfg에 RUN_DATE_ISO 있으면 그 날짜
    #   - 없으면 bry window 중간 날짜
    center_yyyymmdd = pick_run_date(cfg, bry_st + (bry_ed - bry_st)/2)
    center_day = dt.strptime(center_yyyymmdd, "%Y%m%d")
    center_00 = center_day.replace(hour=0, minute=0, second=0, microsecond=0)

    # 최종 window: -6 ~ +30 (3시간 간격)
    t_st = center_00 + timedelta(hours=WINDOW_START_OFFSET_H)
    t_ed = center_00 + timedelta(hours=WINDOW_END_OFFSET_H)

    base_dir = cfg["output"]["base_dir"]
    out_dir = os.path.join(base_dir, "gfs")
    raw_dir = os.path.join(out_dir, "raw")
    thin_dir = os.path.join(out_dir, "thin")
    tmpnc_dir = os.path.join(out_dir, "tmpnc")

    ensure_dir(out_dir)
    ensure_dir(raw_dir)
    ensure_dir(thin_dir)
    ensure_dir(tmpnc_dir)

    out_path = os.path.join(out_dir, f"gfs_korea_{center_yyyymmdd}.nc")

    print("[INFO] config:", CFG_PATH)
    print("[INFO] center day (00Z):", center_00)
    print("[INFO] target window:", t_st, "to", t_ed)
    print("[INFO] dt_hours:", DT_HOURS)
    print("[INFO] output:", out_path)
    print("[INFO] raw:", raw_dir)
    print("[INFO] thin:", thin_dir)

    # build target times
    targets: List[dt] = []
    cur = t_st
    while cur <= t_ed:
        targets.append(cur)
        cur += timedelta(hours=DT_HOURS)

    # 캐시: (run_stamp, fhr) -> raw_path
    raw_cache: Dict[Tuple[str,int], str] = {}

    def get_raw(run_time: dt, fhr: int) -> str:
        stamp = run_time.strftime("%Y%m%d%H")
        key = (stamp, fhr)
        if key in raw_cache and os.path.exists(raw_cache[key]):
            return raw_cache[key]

        out = os.path.join(raw_dir, f"gfs_{stamp}_f{fhr:03d}.grib2")
        if os.path.exists(out):
            raw_cache[key] = out
            return out

        urls = candidate_urls(run_time, fhr)
        print(f"[INFO] download raw run={stamp} f{fhr:03d}")
        print("       urls:", urls[0], " / ", urls[1])

        def _dl():
            return download_any(urls, out, timeout=180)

        used = retry(_dl, tries=DL_TRIES, wait_sec=DL_WAIT_SEC, name=f"download {stamp} f{fhr:03d}")
        if os.path.getsize(out) < 1024 * 200:
            raise RuntimeError(f"downloaded file too small: {out} (url={used})")
        raw_cache[key] = out
        return out

    # per target time: build thin grib -> netcdf -> ds list
    pieces: List[xr.Dataset] = []

    for t in targets:
        run_inst = floor_to_run(t)
        fhr_inst = fhr_from(run_inst, t)
        if fhr_inst not in (0, 3, 6, 9):
            raise RuntimeError(f"instant fhr not in {{0,3,6,9}}: t={t} run={run_inst} fhr={fhr_inst}")

        # cumulative is for interval (t-3, t]
        t_prev = t - timedelta(hours=3)
        run_cum = floor_to_run(t_prev)
        fhr_cum = fhr_from(run_cum, t)
        # 예: t=00, t_prev=21 -> run_cum=18, fhr_cum=6 (3-6 interval ending at 00) 이런 식
        if fhr_cum not in (3, 6, 9):
            raise RuntimeError(f"cumulative fhr not in {{3,6,9}}: t={t} run={run_cum} fhr={fhr_cum}")

        tag = t.strftime("%Y%m%d%H")
        thin_grib = os.path.join(thin_dir, f"gfs_korea_{tag}_thin.grib2")
        tmp_nc    = os.path.join(tmpnc_dir, f"gfs_korea_{tag}.nc")

        print(f"\n[INFO] target={tag}")
        print(f"       instant: run={run_inst.strftime('%Y%m%d%H')} fhr={fhr_inst:03d}")
        print(f"       cum    : run={run_cum.strftime('%Y%m%d%H')} fhr={fhr_cum:03d} (interval {fhr_cum-3}-{fhr_cum})")

        # download raw
        raw_inst = get_raw(run_inst, fhr_inst)
        raw_cum  = get_raw(run_cum,  fhr_cum)

        # build thin (overwrite each time)
        if os.path.exists(thin_grib):
            os.remove(thin_grib)

        inv_inst = wgrib2_lines(raw_inst)
        inv_cum  = wgrib2_lines(raw_cum)

        inst_lines, miss_i = select_instant(inv_inst)
        cum_lines,  miss_c = select_cumulative(inv_cum, fhr_cum)

        if miss_i:
            raise RuntimeError(f"missing instant fields for {tag}: {miss_i}")
        if miss_c:
            raise RuntimeError(f"missing cumulative fields for {tag}: {miss_c}")

        wgrib2_extract(raw_inst, thin_grib, inst_lines)
        wgrib2_append(raw_cum,  thin_grib, cum_lines)

        # thin -> netcdf
        if os.path.exists(tmp_nc):
            os.remove(tmp_nc)
        wgrib2_to_netcdf(thin_grib, tmp_nc)

        ds = xr.open_dataset(tmp_nc)
        ds = force_single_time_and_set(ds, t)
        pieces.append(ds)

    # merge to single
    out = xr.concat(pieces, dim="time").sortby("time")

    tmp = out_path + ".tmp"
    out.to_netcdf(tmp)
    out.close()
    for ds in pieces:
        ds.close()
    os.replace(tmp, out_path)

    print("\n[OK] saved single window file:", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
