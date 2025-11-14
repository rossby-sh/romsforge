
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
avg_surface_ol.py  오버랩(4일) 스킵 + 표층 선택 + 월평균 + 15일 정규화 + 안전 저장

규칙(파일명 끝 4자리/4자리만 사용):
  ..._<cycle:4d>_<chunk:4d>.nc
예: wnp8km_avg_9049_0007.nc

선택 로직(중복 제거):
  - 가장 이른 cycle: 0001~0008 모두 사용
  - 이후 cycle들: 0005~0008만 사용  (기본 overlap=4, chunks_per_cycle=8)

config 예시(config_proc_nft.yaml):
  grdname: /path/to/grid.nc
  data_dir: /path/to/input/
  base_dir: /path/to/output/
  data_prefix_header: wnp8km_avg_     # 비워도 됨("")
  case: 9049
  vars: [zeta, u, v, temp, salt, chlorophyll, NO3, phytoplankton]
  engine: netcdf4                     # 혼재 대응 기본값: netcdf4 (필요시 h5netcdf)
  chunks_per_cycle: 8
  overlap: 4
  validate: true                      # 파일 유효성 검사 on/off
  var_names:
    time_candidates: [ocean_time, time, t]
    lat_candidates:  [lat_rho, lat, y, eta, eta_rho]
    lon_candidates:  [lon_rho, lon, x, xi, xi_rho]
    z_candidates:    [s_rho, z, depth, lev, layer]
    ref_time: ocean_time
"""

# --------- 환경 패치 (NFS/HDF5 안전) ---------
import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# ---------------- 임포트 ----------------
import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import xarray as xr
from dask import compute
from netCDF4 import date2num, Dataset as NC
import pandas as pd

# --- libs/utils.py 로드 ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl  # noqa: E402

# ---------------- [01] 설정 로드 ----------------
cfg = tl.parse_config("./config_proc_nft.yaml")
_ = tl.load_roms_grid(getattr(cfg, "grdname", None))  # 필요 시 확인용

case = getattr(cfg, "case", "case")
data_dir = Path(getattr(cfg, "data_dir", ".")).expanduser().resolve()
base_dir = Path(getattr(cfg, "base_dir", ".")).expanduser().resolve()
outdir = base_dir / str(case)
outdir.mkdir(parents=True, exist_ok=True)

data_prefix_header = str(getattr(cfg, "data_prefix_header", ""))  # 비워도 허용
VARS = list(getattr(cfg, "vars", ["temp", "salt", "zeta"]))
engine = str(getattr(cfg, "engine", "netcdf4")).lower()          # 기본 netcdf4 권장
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

# ---------------- [03] 파일명 파싱 / 오버랩 제거 ----------------
# 접두어 무관: 파일 끝의 4자리/4자리만 파싱
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
        return all_files  # 규칙과 다르면 전체 사용

    cycles = sorted(groups)              # 오름차순
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
    sample = sorted(cnt.items())[:5]
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

# ---------------- [05] 표층 전처리 ----------------
def detect_zdim(da: xr.DataArray):
    # 1) 이름 후보 우선
    for cand in z_candidates:
        if cand in da.dims and da.sizes[cand] > 1:
            return cand
    # 2) 휴리스틱
    not_z = set(time_candidates) | lat_candidates | lon_candidates
    for d in da.dims:
        if d not in not_z and da.sizes[d] > 1:
            return d
    return None  # 2D 변수일 수 있음

def preprocess_surface(ds: xr.Dataset) -> xr.Dataset:
    keep = [v for v in VARS if v in ds.data_vars]
    if not keep:
        return xr.Dataset()
    out = []
    for v in keep:
        da = ds[v]
        zdim = detect_zdim(da)
        out.append(da if zdim is None else da.isel({zdim: -1}, drop=True))
    return xr.merge([x.to_dataset(name=v) for x, v in zip(out, keep)])

# ---------------- [06] 병합 ----------------
ds = xr.open_mfdataset(
    files,
    combine="by_coords",
    preprocess=preprocess_surface,
    parallel=False,             # 안정성 위해 병렬 OFF
    data_vars="minimal",
    coords="minimal",
    compat="override",
    engine=engine,              # 기본 netcdf4 권장
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

# 시간 정렬 + 청크
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
              or "days since 2000-01-01 00:00:00")
time_cal = (src_time.attrs.get("calendar")
            or src_time.encoding.get("calendar")
            or "proleptic_gregorian")

coord_vals = m[time_name].values

if np.issubdtype(m[time_name].dtype, np.datetime64):
    # numpy datetime64 → 15일로 교체
    tvals = pd.to_datetime(coord_vals)
    mid = [pd.Timestamp(year=t.year, month=t.month, day=15, hour=0, minute=0, second=0) for t in tvals]
    m = m.assign_coords({time_name: np.array(mid, dtype="datetime64[ns]")})

    # ★ datetime64일 때는 units/calendar를 attrs/encoding에서 모두 제거
    m[time_name].attrs.pop("units", None)
    m[time_name].attrs.pop("calendar", None)
    m[time_name].encoding.pop("units", None)
    m[time_name].encoding.pop("calendar", None)

else:
    # cftime/수치형 → 15일을 수치로 변환해 대입
    mid_dates = []
    cftime_cls = type(coord_vals[0])
    for t in coord_vals:
        mid_dates.append(cftime_cls(int(t.year), int(t.month), 15, 0, 0, 0))
    mid_nums = np.asarray(date2num(mid_dates, units=time_units, calendar=time_cal), dtype="float64")
    m = m.assign_coords({time_name: (time_name, mid_nums)})

    # ★ 수치형일 땐 attrs 말고 encoding에만 세팅
    m[time_name].attrs.pop("units", None)
    m[time_name].attrs.pop("calendar", None)
    m[time_name].encoding["units"] = time_units
    m[time_name].encoding["calendar"] = time_cal

# 공통 메타(충돌 없는 항목만)
m[time_name].attrs.update({
    "long_name": "time since initialization",
    "field": "time, scalar, series",
    "axis": "T",
})
m[time_name].encoding["_FillValue"] = None
m[time_name].attrs.pop("_FillValue", None)

# 데이터 변수만 인코딩 정리 (좌표는 건드리지 않음)
for v in m.data_vars:
    m[v].encoding.pop("units", None)
    m[v].encoding.pop("calendar", None)

# ---------------- [09] 저장 ----------------
outfile = outdir / f"surface_monthly_avg_{case}.nc"
encoding = {v: {"zlib": True, "complevel": 4, "shuffle": True} for v in m.data_vars}

delayed = m.to_netcdf(outfile, encoding=encoding, engine=engine, compute=False)
compute(delayed)
print(f"[done] → {outfile}")
