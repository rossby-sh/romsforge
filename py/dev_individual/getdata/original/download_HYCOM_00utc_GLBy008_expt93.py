from datetime import datetime as dt, timedelta
import xarray as xr
from netCDF4 import num2date
import numpy as np
from cftime import num2date as cf_num2date

# ====================== 설정 ======================
w_pth = '/home/shjo/DATA/'
d_st = dt(2022, 12, 31)
d_ed = dt(2022, 12, 31)

lat_range = (5, 60)
lon_range = (100, 180)

# HYCOM URL (2023 기준: expt_93.0)
url_ssh = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ssh/2022'
url_ts  = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z/2022'
url_uv   = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z/2022'

# ====================== 데이터 열기 ======================
print('!!! Reading HYCOM expt_93.0 (2024) dataset... !!!')
SSH = xr.open_dataset(url_ssh, decode_times=False)
TEMP = xr.open_dataset(url_ts, decode_times=False)
SALT = xr.open_dataset(url_ts, decode_times=False)
U = xr.open_dataset(url_uv, decode_times=False)
V = xr.open_dataset(url_uv, decode_times=False)

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

        # days since 기반 time 좌표 설정
        time_val = (cur_date - ref_time).total_seconds() / 86400.0
        merged = merged.expand_dims(time=[0])
        merged = merged.assign_coords(time=("time", [time_val]))
        merged['time'].attrs['units'] = 'days since 2000-01-01 00:00:00'
        merged['time'].attrs['calendar'] = 'proleptic_gregorian'

        merged.to_netcdf(
            w_pth + f'HYCOM_{date_str}_00UTC.nc',
            encoding={
                'surf_el': {'zlib': True, 'complevel': 4},
                'water_temp': {'zlib': True, 'complevel': 4},
                'salinity': {'zlib': True, 'complevel': 4},
                'water_u': {'zlib': True, 'complevel': 4},
                'water_v': {'zlib': True, 'complevel': 4},
            },
            unlimited_dims=['time']
        )

        print(f'  Saved HYCOM_{date_str}_00UTC.nc')

    except Exception as e:
        print(f'  Failed on {date_str}: {e}')

    cur_date += timedelta(days=1)

