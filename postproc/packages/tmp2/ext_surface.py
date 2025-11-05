
import os
from pathlib import Path
import xarray as xr
from dask import compute
import sys
from netCDF4 import Dataset
import numpy as np

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
print(f"[info] 대상 변수(표층만 저장): {VARS}")

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

# --- [03] 표층만 유지해서 열기 (시간축/원래 units 그대로) ---
ds = xr.open_mfdataset(
    files,
    combine="by_coords",
    preprocess=preprocess_surface,
    parallel=False,
    data_vars="minimal",
    coords="minimal",
    compat="override",
    engine="netcdf4",   # 읽기는 netCDF4
)

# 시간좌표 결정
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
print(f"[info] using time coord: {time_name}")

# 안정성: time 청크 1
ds = ds.chunk({time_name: 1})

# 표층 Dataset 그대로 사용
m = ds

# --- [03-1] ncview/ROMS 호환 메타(시간축) 보강 ---
# 원본 units/calendar는 그대로 두고, field/axis만 추가
m[time_name].attrs.setdefault("field", "time, scalar, series")
m[time_name].attrs.setdefault("axis", "T")

# 좌표의 _FillValue 제거 (좌표엔 불필요)
m[time_name].encoding["_FillValue"] = None
m[time_name].attrs.pop("_FillValue", None)

# (중요) 인코딩 정리는 데이터 변수만 (좌표는 건드리지 않음)
for v in m.data_vars:
    m[v].encoding.pop("units", None)
    m[v].encoding.pop("calendar", None)

# 청크 정리 후 메모리에 로드(파일 핸들 의존 제거)
m = m.unify_chunks().load()

# --- [04] 저장 (netCDF4, unlimited time 차원) ---
outfile = outdir / f"surface_timeseries_surfaceonly_{case}.nc"
encoding = {v: {"zlib": True, "complevel": 4, "shuffle": True} for v in m.data_vars}
unlim = {time_name}

t = m.to_netcdf(
    outfile,
    encoding=encoding,
    engine="netcdf4",
    format="NETCDF4_CLASSIC",
    unlimited_dims=unlim,
    compute=False,
)
compute(t)
print(f"done. → {outfile}")
