import os
from pathlib import Path
import xarray as xr
from dask import compute

pth = "/data/shjo/applications/fennel_version_test/storage/roms_fennel_4_0_org/"
files = sorted([str(Path(pth)/f) for f in os.listdir(pth) if f.startswith("4_0_org_avg")])

print(files)
VARS = ["temp", "salt", "chlorophyll", "NO3"]

# 수직차원 자동 감지 (time/lat/lon 제외)
def detect_zdim(var):
    # 자주 쓰는 좌표 이름들
    not_z = {"ocean_time","t","lat","latitude","y","eta","eta_rho",
             "lon","longitude","x","xi","xi_rho"}
    for d in var.dims:
        if d not in not_z and var.sizes[d] > 1:
            return d
    # 못 찾으면 흔한 이름 중 존재하는 것 선택
    for cand in ("s_rho","z","depth","lev","layer"):
        if cand in var.dims:
            return cand
    raise ValueError(f"수직 차원을 찾지 못했어: dims={var.dims}")

# 필요한 변수만 남기고, 표층(z=-1)만 잘라서 반환
def preprocess_surface(ds):
    keep = [v for v in VARS if v in ds.data_vars]
    ds = ds[keep]
    # 대표 변수로 zdim 판단 (temp 우선)
    base = ds[keep[0]]
    zdim = detect_zdim(base)
    # 모든 변수 표층만 선택 (time은 그대로 유지)
    sel = {zdim: -1}
    ds = ds.isel(sel, drop=True)  # drop=True로 zdim 제거 → 3D(time, lat, lon)
    return ds

# 한 번만 열고, 파일당 time=1 구조 가정 → chunks={'time':1} 권장
ds = xr.open_mfdataset(
    files,
    combine="by_coords",
    preprocess=preprocess_surface,
    parallel=True,
    chunks={"ocean_time": 1},
    data_vars="minimal",
    coords="minimal",
    compat="override",
    # engine="h5netcdf",  # 필요시 엔진 고정
)

# 월평균(월초 기준). 표층만 남겨놨으므로 3D: (time, lat, lon)
print("monthly mean...")
temp_m = ds["temp"].resample(ocean_time="MS").mean()
salt_m = ds["salt"].resample(ocean_time="MS").mean()
chl_m  = ds["chlorophyll"].resample(ocean_time="MS").mean()
no3_m  = ds["NO3"].resample(ocean_time="MS").mean()

# 지연(Defer) 저장 태스크를 만들고 한 번에 compute → 오버헤드 절감
t1 = temp_m.to_netcdf("./temp_monthly.nc", compute=False)
t2 = salt_m.to_netcdf("./salt_monthly.nc", compute=False)
t3 = chl_m.to_netcdf("./chl_monthly.nc", compute=False)
t4 = no3_m.to_netcdf("./no3_monthly.nc", compute=False)
compute(t1, t2, t3, t4)
print("done.")
