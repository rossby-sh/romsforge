
import os
from pathlib import Path
import xarray as xr
from dask import compute
import sys
from netCDF4 import Dataset, date2num
import numpy as np
import pandas as pd
import datetime as dt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..', 'libs')))
import utils as tl

# --- [01] Load configuration ---
cfg = tl.parse_config("./config_proc.yaml")
grd = tl.load_roms_grid(cfg.grdname)

case = getattr(cfg, "case", "case")
data_dir = Path(getattr(cfg, "data_dir", ".")).expanduser().resolve()
data_prefix_header = getattr(cfg, "data_prefix_header", "")
base_dir = Path(getattr(cfg, "base_dir", ".")).expanduser().resolve()
VARS = getattr(cfg, "vars", ["temp", "salt", "chlorophyll", "NO3"])

# var_names 블록 파싱
vn = getattr(cfg, "var_names", None) or {}
time_candidates = list(getattr(vn, "time_candidates", ["ocean_time", "time", "t"]))
lat_candidates  = set(getattr(vn, "lat_candidates",  ["lat", "latitude", "y", "eta", "eta_rho"]))
lon_candidates  = set(getattr(vn, "lon_candidates",  ["lon", "longitude", "x", "xi", "xi_rho"]))
z_candidates    = list(getattr(vn, "z_candidates",   ["s_rho", "z", "depth", "lev", "layer"]))
ref_time_name   = getattr(vn, "ref_time", None)  # 우선 사용할 시간 좌표 이름

# 출력 디렉토리
outdir = base_dir / case
outdir.mkdir(parents=True, exist_ok=True)

# 입력 파일 목록
files = sorted([
    str(p) for p in data_dir.iterdir()
    if p.is_file() and p.name.startswith(data_prefix_header)
       and p.suffix.lower() in {".nc", ".nc4", ".cdf"}
])
if not files:
    raise FileNotFoundError(f"입력 파일 없음: {data_dir}, prefix={data_prefix_header}")

print(f"[info] files={len(files)}개, 예: {files[:3]}")
print(f"[info] 평균 변수: {VARS}")

# --- [02] 표층 선택용 함수 ---
def detect_zdim(var):
    not_z = set(time_candidates) | lat_candidates | lon_candidates
    for d in var.dims:
        if d not in not_z and var.sizes[d] > 1:
            return d
    for cand in z_candidates:
        if cand in var.dims:
            return cand
    raise ValueError(f"수직 차원을 찾지 못했어: dims={var.dims}")

def preprocess_surface(ds):
    keep = [v for v in VARS if v in ds.data_vars]
    if not keep:
        return xr.Dataset()
    ds = ds[keep]
    zdim = detect_zdim(ds[keep[0]])
    return ds.isel({zdim: -1}, drop=True)  # 표층만

# --- [03] 병합 후 월평균 ---
ds = xr.open_mfdataset(
    files,
    combine="by_coords",
    preprocess=preprocess_surface,
    parallel=True,
    data_vars="minimal",
    coords="minimal",
    compat="override",
    engine="netcdf4",   # 읽기도 h5netcdf로 고정
)

# 시간좌표 결정 (ref_time_name 우선, 없으면 후보 순회)
time_name = None
if ref_time_name and ref_time_name in ds.coords:
    time_name = ref_time_name
else:
    for cand in time_candidates:
        if cand in ds.coords:
            time_name = cand
            break
if time_name is None:
    raise KeyError(f"시간 좌표 후보({time_candidates}) 중 찾은 게 없어.")
print(f"[info] monthly mean using time: {time_name}")

# 시간축 청크 적용(안정성)
ds = ds.chunk({time_name: 1})

# 평균 계산
merged = []
for v in VARS:
    if v not in ds:
        print(f"[warn] {v} 없음 → 스킵")
        continue
    avg = ds[v].resample({time_name: "MS"}).mean()  # 월초 기준
    merged.append(avg.to_dataset(name=v))
if not merged:
    raise RuntimeError("평균 낼 변수가 하나도 없어.")
m = xr.merge(merged)

# --- [03-1] 월평균 시간좌표를 '그 달 15일 00:00'로 맞춤 ---
src_time = ds[time_name]
time_units = (src_time.attrs.get("units")
              or src_time.encoding.get("units")
              or "days since 1970-01-01 00:00:00")
time_cal = (src_time.attrs.get("calendar")
            or src_time.encoding.get("calendar")
            or "standard")

coord_vals = m[time_name].values
mid_dates = []
if np.issubdtype(m[time_name].dtype, np.datetime64):
    for t in pd.to_datetime(coord_vals):
        mid_dates.append(dt.datetime(int(t.year), int(t.month), 15, 0, 0, 0))
else:
    cftime_cls = type(coord_vals[0])
    for t in coord_vals:
        mid_dates.append(cftime_cls(int(t.year), int(t.month), 15, 0, 0, 0))

mid_nums = np.asarray(date2num(mid_dates, units=time_units, calendar=time_cal), dtype="float64")
# (time_name, 값) 형태로 차원-좌표 변수로 명시
m = m.assign_coords({time_name: (time_name, mid_nums)})

# 시간좌표 attrs 유지(encoding엔 넣지 말 것)
m[time_name].attrs["units"] = time_units
m[time_name].attrs["calendar"] = time_cal

# 모든 변수/좌표의 encoding에서 units/calendar 제거
for name in list(m.variables):
    m[name].encoding.pop("units", None)
    m[name].encoding.pop("calendar", None)

# 청크 정리 후 메모리에 로드(파일 핸들 의존 제거)
m = m.unify_chunks()
m = m.load()

m[time_name].attrs.update({
    "long_name": "time since initialization",
    "calendar": "gregorian",
    "field": "time, scalar, series",
    "axis": "T",
})
m[time_name].encoding["_FillValue"] = None
m[time_name].attrs.pop("_FillValue", None)

# --- [04] 저장 ---
outfile = outdir / f"surface_monthly_avg_{case}.nc"
# 데이터 변수만 압축 인코딩 (좌표 X)
encoding = {v: {"zlib": True, "complevel": 4, "shuffle": True} for v in m.data_vars}

t = m.to_netcdf(outfile, encoding=encoding, engine="netcdf4", compute=False)
compute(t)
print(f"done. → {outfile}")
