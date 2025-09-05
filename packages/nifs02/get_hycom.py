from datetime import datetime as dt, timedelta
import xarray as xr
from netCDF4 import num2date
import numpy as np
from cftime import num2date as cf_num2date
import yaml
import os
# ====================== 설정 ======================

with open("config_all.yaml") as f:
    cfg = yaml.safe_load(f)

# config에서 불러오기
w_pth = cfg["output"]["base_dir"]+"hycom/"
os.makedirs(w_pth, exist_ok=True)
d_st = dt.strptime(cfg["time"]["start"], "%Y-%m-%d")
d_ed = dt.strptime(cfg["time"]["end"],   "%Y-%m-%d")
lat_range = tuple(cfg["region"]["lat"])
lon_range = tuple(cfg["region"]["lon"])

print(d_st)
print(d_ed)
print(lon_range)
print(w_pth)

print("=== Download HYCOM ===")
print("--- Save to directory: "+w_pth+" ---")

# HYCOM URL (2024+)
url_ssh  = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/ssh'
url_temp = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/t3z'
url_salt = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/s3z'
url_u    = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/u3z'
url_v    = 'https://tds.hycom.org/thredds/dodsC/ESPC-D-V02/v3z'

# ====================== 데이터 열기 ======================
print('!!! Reading HYCOM dataset... !!!')
SSH = xr.open_dataset(url_ssh, decode_times=False)
TEMP = xr.open_dataset(url_temp, decode_times=False)
SALT = xr.open_dataset(url_salt, decode_times=False)
U = xr.open_dataset(url_u, decode_times=False)
V = xr.open_dataset(url_v, decode_times=False)

units = 'hours since 2000-01-01 00:00:00'
ssh_times = np.array(num2date(SSH.time.values, units))
var_times = np.array(num2date(TEMP.time.values, units))
ssh_dates = np.array([d.strftime('%Y-%m-%d') for d in ssh_times])
var_dates = np.array([d.strftime('%Y-%m-%d') for d in var_times])

# ====================== 루프 시작 ======================
cur_date = d_st
ref_time = dt(2000, 1, 1)

while cur_date <= d_ed:
    date_str = cur_date.strftime('%Y%m%d')
    print(f'\n--- Processing {date_str} ---')

    cur_str = cur_date.strftime('%Y-%m-%d')
    ssh_idx = np.where((ssh_dates == cur_str) & (np.array([d.hour for d in ssh_times]) == 0))[0]
    var_idx = np.where((var_dates == cur_str) & (np.array([d.hour for d in var_times]) == 0))[0]

    if len(ssh_idx) == 0 or len(var_idx) == 0:
        print(f'  Skipping {date_str}: no 00 UTC data')
        cur_date += timedelta(days=1)
        continue

    try:
        subset_ssh = SSH['surf_el'].isel(time=ssh_idx[0]).sel(lat=slice(*lat_range), lon=slice(*lon_range))
        subset_temp = TEMP['water_temp'].isel(time=var_idx[0]).sel(lat=slice(*lat_range), lon=slice(*lon_range))
        subset_salt = SALT['salinity'].isel(time=var_idx[0]).sel(lat=slice(*lat_range), lon=slice(*lon_range))
        subset_u = U['water_u'].isel(time=var_idx[0]).sel(lat=slice(*lat_range), lon=slice(*lon_range))
        subset_v = V['water_v'].isel(time=var_idx[0]).sel(lat=slice(*lat_range), lon=slice(*lon_range))

        merged = xr.merge([subset_ssh, subset_temp, subset_salt, subset_u, subset_v])

        # days since 기반 time 좌표 지정
        time_val = (cur_date - ref_time).total_seconds() / 86400.0
        merged = merged.expand_dims(time=[0])
        merged = merged.assign_coords(time=("time", [time_val]))
        merged['time'].attrs['units'] = 'days since 2000-01-01 00:00:00'
        merged['time'].attrs['calendar'] = 'proleptic_gregorian'

        merged.to_netcdf(
            w_pth + f'hycom_{date_str}_00utc.nc',
            encoding={
                'surf_el': {'zlib': True, 'complevel': 4},
                'water_temp': {'zlib': True, 'complevel': 4},
                'salinity': {'zlib': True, 'complevel': 4},
                'water_u': {'zlib': True, 'complevel': 4},
                'water_v': {'zlib': True, 'complevel': 4},
            },
            unlimited_dims=['time']
        )

        print(f'  Saved {date_str}_00utc.nc')

    except Exception as e:
        print(f'  Failed on {date_str}: {e}')

    cur_date += timedelta(days=1)



