# -*- coding: utf-8 -*-
"""
Created on Thu May 15 13:31:16 2025

@author: ust21
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata

# 파일 열기
obs_file = 'D:/shjo/ROMS_OUTS/E_test/roms_obs_phyt_30km_N36_250501_250531_03E_.nc'
obs_file = 'D:/shjo/roms_obs_phyt_30km_N36_250501_250531_30SE_.nc'

ds = nc.Dataset(obs_file)

# 변수 추출
lon = ds.variables['obs_lon'][:]
lat = ds.variables['obs_lat'][:]
val = ds.variables['obs_value'][:]
err = ds.variables['obs_error'][:]
time = ds.variables['obs_time'][:]
type = ds.variables['obs_type'][:] if 'obs_type' in ds.variables else None

# 단위 변환 (옵션): 시간 단위를 문자열로 확인
time_unit = ds.variables['obs_time'].units
time_unit = "days since 2000-1-1"
calendar = ds.variables['obs_time'].calendar if 'calendar' in ds.variables['obs_time'].ncattrs() else 'standard'
time_dates = nc.num2date(time, time_unit, calendar=calendar)

# 예: SST (type == 1)만 필터
if type is not None:
    sst_idx = np.where(type == 9)
    lon, lat, val, err, time_dates = lon[sst_idx], lat[sst_idx], val[sst_idx], err[sst_idx], time_dates[sst_idx]

# 고정 그리드 생성 (해역 범위에 따라 조정)
lon_min, lon_max = lon.min(), lon.max()
lat_min, lat_max = lat.min(), lat.max()
grid_x, grid_y = np.meshgrid(
    np.linspace(lon_min, lon_max, 200),
    np.linspace(lat_min, lat_max, 200)
)


# 고유 시간 구간으로 반복
unique_times = np.unique(time_dates)



for t in unique_times:
    idx = np.where(time_dates == t)

    grid_val = griddata(
    (lon[idx], lat[idx]), err[idx],
    (grid_x, grid_y), method='linear'  # or 'cubic', 'nearest'
    )

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_title(f'Phyt Observation Error - {t.strftime("%Y-%m-%d %H:%M:%S")}', fontsize=14)
    ax.set_extent([105, 165, 15, 50], crs=ccrs.PlateCarree())

    # 해안선 등 추가
    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='#ffffff')
    ax.gridlines(draw_labels=True, linestyle='--')

    # 관측 산점도
    val_max=np.max(err[idx])
    idV=np.argmax(err[idx])
    
    lonV=lon[idV]
    latV=lat[idV]

    ax.text(
    0.01, 0.99,  # (x, y): 좌측 상단 기준 (axes fraction)
    f'Max: {val_max:.2f} mmol/m^3\n'+f'Max lon: {lonV:.2f}\n'+f'Max lat: {latV:.2f}' ,
    transform=ax.transAxes,  # 축 좌표계 기준
    fontsize=12,
    fontweight='bold',
    va='top',
    ha='left',
    bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')  # 배경 박스 옵션
)
    sc = ax.scatter(lon[idx], lat[idx], c=err[idx], cmap='jet', s=10,vmin=0, vmax=0.05,
                    edgecolor=None, transform=ccrs.PlateCarree())

    # 컬러바
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.06, shrink=.8)
    cbar.set_label('Phyt Observation error (mmol m-3)', fontsize=12)

    plt.tight_layout()
    # plt.savefig(f'D:/shjo/ROMS_OUTS/E_test/Case_03E/figs/obs_error/phyto_err_{t.strftime("%Y%m%d")}', dpi=200)

    plt.show()
    
    
    
    
    
    
    