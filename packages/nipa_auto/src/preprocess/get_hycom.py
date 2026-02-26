
from datetime import datetime as dt, timedelta
import xarray as xr
from netCDF4 import num2date
import numpy as np
import yaml
import os

# ====================== 설정 ======================
with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

# config에서 불러오기
w_pth = os.path.join(cfg["output"]["base_dir"], "hycom")
os.makedirs(w_pth, exist_ok=True)

# 수정
t_st = cfg["bry_start_date"]
t_ed = cfg["bry_end_date"]

if isinstance(t_st, str):
    t_st = dt.strptime(t_st, "%Y-%m-%d %H:%M:%S")
if isinstance(t_ed, str):
    t_ed = dt.strptime(t_ed, "%Y-%m-%d %H:%M:%S")



lat_range = tuple(cfg["region"]["lat"])
lon_range = tuple(cfg["region"]["lon"])

# lat/lon slice 방향 안전하게
lat1, lat2 = lat_range
lon1, lon2 = lon_range
lat_slice = slice(min(lat1, lat2), max(lat1, lat2))
lon_slice = slice(min(lon1, lon2), max(lon1, lon2))

print(t_st, t_ed)
print(lat_range, lon_range)
print(w_pth)

print("=== Download HYCOM ===")
print("--- Save to directory:", w_pth, "---")

# HYCOM URL (2024+)
url_ssh  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/ssh'
url_temp = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/t3z'
url_salt = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/s3z'
url_u    = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/u3z'
url_v    = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/v3z'

# ====================== 데이터 열기 ======================
print('!!! Reading HYCOM dataset... !!!')
SSH  = xr.open_dataset(url_ssh,  decode_times=False)
TEMP = xr.open_dataset(url_temp, decode_times=False)
SALT = xr.open_dataset(url_salt, decode_times=False)
U    = xr.open_dataset(url_u,    decode_times=False)
V    = xr.open_dataset(url_v,    decode_times=False)

units = 'hours since 2000-01-01 00:00:00'
ref_time = dt(2000, 1, 1)

# time decode (python datetime array)
ssh_times = np.array(num2date(SSH.time.values, units))
var_times = np.array(num2date(TEMP.time.values, units))  # TEMP/SALT/U/V가 같은 time이라고 가정

# 빠른 인덱싱용: datetime -> index
ssh_map = {t: i for i, t in enumerate(ssh_times)}
var_map = {t: i for i, t in enumerate(var_times)}

# ====================== 3시간 루프 ======================
cur = t_st
step = timedelta(hours=3)

while cur <= t_ed:
    # 저장 이름용(UTC 기준이라고 전제)
    stamp = cur.strftime("%Y%m%d%H")
    print(f"\n--- Processing {cur.strftime('%F %T')} -> {stamp} ---")

    # 정확히 일치하는 time index 찾기
    ssh_i = ssh_map.get(cur, None)
    var_i = var_map.get(cur, None)

    if ssh_i is None:
        print(f"  Skip {stamp}: no SSH data at this time")
        cur += step
        continue
    if var_i is None:
        print(f"  Skip {stamp}: no 3D var data at this time")
        cur += step
        continue

    try:
        subset_ssh  = SSH['surf_el'].isel(time=ssh_i).sel(lat=lat_slice, lon=lon_slice)
        subset_temp = TEMP['water_temp'].isel(time=var_i).sel(lat=lat_slice, lon=lon_slice)
        subset_salt = SALT['salinity'].isel(time=var_i).sel(lat=lat_slice, lon=lon_slice)
        subset_u    = U['water_u'].isel(time=var_i).sel(lat=lat_slice, lon=lon_slice)
        subset_v    = V['water_v'].isel(time=var_i).sel(lat=lat_slice, lon=lon_slice)

        merged = xr.merge([subset_ssh, subset_temp, subset_salt, subset_u, subset_v])

        # time 좌표: days since 2000-01-01 (시간 포함, fraction)
        time_val = (cur - ref_time).total_seconds() / 86400.0
        merged = merged.expand_dims(time=[0])
        merged = merged.assign_coords(time=("time", [time_val]))
        merged['time'].attrs['units'] = 'days since 2000-01-01 00:00:00'
        merged['time'].attrs['calendar'] = 'proleptic_gregorian'

        out_fn = os.path.join(w_pth, f"hycom_korea_{stamp}.nc")
        merged.to_netcdf(
            out_fn,
            encoding={
                'surf_el': {'zlib': True, 'complevel': 4},
                'water_temp': {'zlib': True, 'complevel': 4},
                'salinity': {'zlib': True, 'complevel': 4},
                'water_u': {'zlib': True, 'complevel': 4},
                'water_v': {'zlib': True, 'complevel': 4},
            },
            unlimited_dims=['time']
        )

        print(f"  Saved {out_fn}")

    except Exception as e:
        print(f"  Failed on {stamp}: {e}")

    cur += step
