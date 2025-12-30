
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
avg_surface_ol.py  오버랩(4일) 스킵 + 전층 유지(3D) + 월평균 + 15일 정규화 (netCDF4 안정 저장)
- 파일명 끝 '<cycle:4d>_<chunk:4d>.nc'만 파싱 (접두어 변화 무관)
- 첫 사이클: 0001~0008, 이후 사이클: 0005~0008만 사용
- 엔진: netcdf4 (열기/저장)
- NFS/HDF5 이슈 회피: HDF5_USE_FILE_LOCKING=FALSE, open_mfdataset 병렬 OFF
- ocean_time이 datetime64여도 encoding에 units/calendar를 명시해 절대 기준 숫자로 저장
"""

# --------- NFS/HDF5 안전 패치 ---------
import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# ---------------- 임포트 ----------------
import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import xarray as xr
import pandas as pd
import dask
from netCDF4 import date2num, Dataset as NC

# --- libs/utils.py 로드 ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl  # noqa: E402


# ---------------- [01] 설정 로드 ----------------
cfg = tl.parse_config("./config_proc_nft.yaml")
_ = tl.load_roms_grid(getattr(cfg, "grdname", None))  # 필요 시 확인용(실제 사용 안 해도 됨)

case = getattr(cfg, "case", "case")
data_dir = Path(getattr(cfg, "data_dir", ".")).expanduser().resolve()
base_dir = Path(getattr(cfg, "base_dir", ".")).expanduser().resolve()
outdir = base_dir / str(case)
outdir.mkdir(parents=True, exist_ok=True)

data_prefix_header = str(getattr(cfg, "data_prefix_header", ""))  # 비워도 허용
VARS = list(getattr(cfg, "vars", ["temp", "salt", "zeta"]))
engine = "netcdf4"                                 # 고정
chunks_per_cycle = int(getattr(cfg, "chunks_per_cycle", 8))
overlap = int(getattr(cfg, "overlap", 4))
validate = bool(getattr(cfg, "validate", True))

vn = getattr(cfg, "var_names", None) or {}
time_candidates = list(getattr(vn, "time_candidates", ["ocean_time", "time", "t"]))
lat_candidates  = set(getattr(vn, "lat_candidates",  ["lat", "latitude", "y", "eta", "eta_rho"]))
lon_candidates  = set(getattr(vn, "lon_candidates",  ["lon", "longitude", "x", "xi", "xi_rho"]))
z_candidates    = list(getattr(vn, "z_candidates",   ["s_rho", "z", "depth", "lev", "layer"]))
ref_time_name   = getattr(vn, "ref_time", None)


# ---------------- [02] 입력 파일 스캔 ----------------
cands = []
for p in data_dir.iterdir():
    if not p.is_file():
        continue
    if p.suffix.lower() not in {".nc", ".nc4", ".cdf"}:
        continue
    if data_prefix_header and not p.name.startswith(data_prefix_header):
        continue
    cands.append(str(p))
files_all = sorted(cands)
if not files_all:
    raise FileNotFoundError(f"입력 파일 없음: {data_dir}, prefix={data_prefix_header!r}")

print(f"[info] files={len(files_all)}개 (예: {files_all[:3]})")
print(f"[info] 평균 변수: {VARS}")


# ---------------- [03] 오버랩 제거 ----------------
NAME_RE = re.compile(r"(?P<cycle>\d{4})_(?P<chunk>\d{4})\.nc$")

def select_nonoverlap(all_files, chunks_per_cycle=8, overlap=4):
    groups = {}  # cycle -> [(chunk, path)]
    for f in all_files:
        m = NAME_RE.search(Path(f).name)
        if not m:
            continue
        cyc = int(m.group("cycle")); chk = int(m.group("chunk"))
        groups.setdefault(cyc, []).append((chk, f))
    if not groups:
        return all_files
    cycles = sorted(groups)
    first_cycle = cycles[0]
    chosen = []
    for cyc in cycles:
        s = 1 if cyc == first_cycle else overlap + 1
        e = chunks_per_cycle
        for chk, f in sorted(groups[cyc]):
            if s <= chk <= e:
                chosen.append(f)
    return chosen

files = select_nonoverlap(files_all, chunks_per_cycle=chunks_per_cycle, overlap=overlap)
print(f"[info] 비중복 선택 후 files={len(files)}개")

# 디버그: 사이클별 선택 개수 요약
cnt = defaultdict(int)
for f in files:
    m = NAME_RE.search(Path(f).name)
    if m:
        cnt[int(m.group("cycle"))] += 1
if cnt:
    sample = sorted(cnt.items())[:10]
    mean_sel = np.mean([v for _, v in cnt.items()])
    print(f"[check] 사이클별 선택 개수(앞5): {sample}  (평균≈{mean_sel:.2f})")


# ---------------- [04] 유효 파일 검증(옵션) ----------------
def filter_valid_netcdf(paths):
    if not validate:
        return paths
    ok, bad = [], []
    for f in paths:
        try:
            with NC(f, "r"):
                pass
            ok.append(f)
        except Exception as e:
            bad.append((f, str(e)))
    if bad:
        print(f"[warn] 손상/비호환 파일 {len(bad)}개 제외 (예시 5개):")
        for f, msg in bad[:5]:
            print("  -", f, "=>", msg)
    return ok

files = filter_valid_netcdf(files)
print(f"[info] 검증 통과 files={len(files)}개")
if not files:
    raise RuntimeError("열 수 있는 파일이 없음.")


# ---------------- [05] 전층 전처리 ----------------
def preprocess_surface(ds: xr.Dataset) -> xr.Dataset:
    """
    VARS에 지정된 변수만 남기고, 수직 차원은 그대로 유지.
    이후 단계에서 시간축을 따라 월평균만 수행.
    """
    keep = [v for v in VARS if v in ds.data_vars]
    if not keep:
        return xr.Dataset()
    return ds[keep]


# ---------------- [06] 병합 (netCDF4, 병렬 OFF) ----------------
ds = xr.open_mfdataset(
    files,
    combine="by_coords",
    preprocess=preprocess_surface,
    parallel=False,             # 중요: libnetcdf/HDF5 충돌 방지
    data_vars="minimal",
    coords="minimal",
    compat="override",
    engine=engine,              # netcdf4
)

# 시간좌표 선택
time_name = None
if ref_time_name and ref_time_name in ds.coords:
    time_name = ref_time_name
else:
    for cand in time_candidates:
        if cand in ds.coords:
            time_name = cand
            break
if time_name is None:
    raise KeyError(f"시간 좌표 후보({time_candidates}) 중 찾지 못함.")
print(f"[info] monthly mean using time: {time_name}")

ds = ds.sortby(time_name)
ds = ds.chunk({time_name: 1})


# ---------------- [07] 월평균 ----------------
merged = []
for v in VARS:
    if v not in ds:
        print(f"[warn] {v} 없음 → 스킵")
        continue
    avg = ds[v].resample({time_name: "MS"}, label="left", closed="left").mean()
    merged.append(avg.to_dataset(name=v))
if not merged:
    raise RuntimeError("평균 낼 변수가 하나도 없어.")
m = xr.merge(merged)


# ---------------- [08] 시간좌표 15일 00:00로 정규화 ----------------
src_time = ds[time_name]
time_units = (src_time.attrs.get("units")
              or src_time.encoding.get("units")
              or "days since 2000-01-01 00:00:00")   # 원하는 기준으로 고정
time_cal = (src_time.attrs.get("calendar")
            or src_time.encoding.get("calendar")
            or "proleptic_gregorian")

coord_vals = m[time_name].values

if np.issubdtype(m[time_name].dtype, np.datetime64):
    # datetime64 → 15일로 교체
    tvals = pd.to_datetime(coord_vals)
    mid = [pd.Timestamp(year=t.year, month=t.month, day=15, hour=0, minute=0, second=0) for t in tvals]
    m = m.assign_coords({time_name: np.array(mid, dtype="datetime64[ns]")})

    # attrs 충돌 제거
    m[time_name].attrs.pop("units", None)
    m[time_name].attrs.pop("calendar", None)
    m[time_name].attrs.pop("_FillValue", None)

    # ★ datetime64라도 encoding에 units/calendar를 명시 → 절대 기준 숫자로 저장됨
    m[time_name].encoding["units"] = time_units
    m[time_name].encoding["calendar"] = time_cal
    m[time_name].encoding["_FillValue"] = None

else:
    # cftime/수치형 → 15일을 수치로 변환
    mid_dates = []
    cftime_cls = type(coord_vals[0])
    for t in coord_vals:
        mid_dates.append(cftime_cls(int(t.year), int(t.month), 15, 0, 0, 0))
    mid_nums = np.asarray(date2num(mid_dates, units=time_units, calendar=time_cal), dtype="float64")
    m = m.assign_coords({time_name: (time_name, mid_nums)})

    # attrs 말고 encoding에만 세팅
    m[time_name].attrs.pop("units", None)
    m[time_name].attrs.pop("calendar", None)
    m[time_name].encoding["units"] = time_units
    m[time_name].encoding["calendar"] = time_cal
    m[time_name].encoding["_FillValue"] = None

# 공통 메타(정보용만 남김)
m[time_name].attrs.update({
    "long_name": "time since initialization",
    "field": "time, scalar, series",
    "axis": "T",
})

# 데이터 변수만 인코딩 정리
for v in m.data_vars:
    m[v].encoding.pop("units", None)
    m[v].encoding.pop("calendar", None)


# ---------------- [09] 저장 (netCDF4, single-threaded) ----------------
outfile = outdir / f"surface_monthly_avg_{case}.nc"
encoding = {v: {"zlib": True, "complevel": 4, "shuffle": True} for v in m.data_vars}

# 단일 스레드로 동기 저장 (교착/세그폴트 회피)
with dask.config.set(scheduler="single-threaded"):
    m.to_netcdf(outfile, encoding=encoding, engine=engine, format="NETCDF4")

print(f"[done] → {outfile}")
