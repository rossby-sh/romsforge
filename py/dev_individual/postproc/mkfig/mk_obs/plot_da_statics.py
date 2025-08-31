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
# obs_file = 'D:/shjo/roms_obs_phyt_30km_N36_250501_250531_30SE_.nc'
obs_file = 'D:/shjo/roms_obs_phyt_30km_N36_250630_250801_.nc'
# obs_file = 'D:/shjo/roms_obs_phyt_jul02_anal.nc'

ds = nc.Dataset(obs_file)

# 변수 추출
lon = ds.variables['obs_lon'][:]
lat = ds.variables['obs_lat'][:]
val = ds.variables['obs_value'][:]
err = ds.variables['obs_error'][:]
xgrid = ds.variables['obs_Xgrid'][:] 
ygrid = ds.variables['obs_Ygrid'][:] 

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
    x, y, val, err, time_dates = xgrid[sst_idx], ygrid[sst_idx], val[sst_idx], err[sst_idx], time_dates[sst_idx]
    lon, lat = lon[sst_idx], lat[sst_idx]

# ncG=nc.Dataset("D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc")
# lon_rho,lat_rho = ncG['lon_rho'][:].data, ncG["lat_rho"][:].data
# val[val>=0.02**2]=np.nan

# lon_rho[int(x.data[0]),int(y.data[0])]

# 고정 그리드 생성 (해역 범위에 따라 조정)
lon_min, lon_max = lon.min(), lon.max()
lat_min, lat_max = lat.min(), lat.max()
grid_x, grid_y = np.meshgrid(
    np.linspace(lon_min, lon_max, 200),
    np.linspace(lat_min, lat_max, 200)
)


# 고유 시간 구간으로 반복
unique_times = np.unique(time_dates)



import datetime as dt

a12=[]; data_len=[]
lower=[]

fig, ax = plt.subplots(figsize=(10, 5))
ax.set_xlim(1, 31)
ax.set_xticks(range(1, 32))  
plt.xlim([0.5,31.5])
n=0
for t in unique_times:
    idx = np.where( time_dates == t )
    N=np.ones(len(idx[0]))*n

    a12.append(np.where( (val.data[idx[0]] < err.data[idx[0]]) & (err.data[idx[0]] >0.02) )[0].shape[0])
    lower.append(np.where( val.data[idx[0]] < 0.02)[0].shape[0])

    
    # plt.plot(np.ones(len(idx[0])),val[idx],'r')
    plt.plot(N,err[idx],'b',marker='o',linewidth=0,markersize=3)
    plt.plot(N,val[idx],'r',marker='o',linewidth=0,markersize=3)
    n+=1
    # a12.append(np.sum(val[idx]>=11))
    data_len.append(len(idx[0]))
    
ax.legend(['obs_error','obs_values'],fontsize=16)
ax.set_xticks(np.arange(1,32,2))
ax.set_xticklabels([str(i) for i in range(1, 32,2)],fontsize=16)
ax.set_yticks(np.arange(0,250,50))
ax.set_yticklabels([str(i) for i in range(0, 250,50)],fontsize=16)

# [str(ii)+'/'+str(jj) for ii,jj in zip(a12,data_len)]

import pandas as pd

AA=pd.DataFrame({},columns=['upper','lower','normal','total','percent'])

AA['upper'] = a12
AA['lower'] = lower
AA['total'] = data_len
AA['normal']=AA['total'] - AA['lower'] - AA['upper']
AA['percent1'] = (AA['upper']+AA['lower'])/AA['total']*100
AA['percent'] = AA['percent1'].map(lambda x: f"{x:.2f}")

AA.to_csv("D:/shjo/jul_obs.csv")






fig, ax = plt.subplots(figsize=(10, 5))
ax.set_xlim(1, 31)
ax.set_xticks(range(1, 32))  
plt.xlim([0.5,31.5])
plt.ylim([0,30])
n=0

ax.plot(AA['percent1'].values,marker='o',color='C3')
ax.axhline(y=10,color='k',alpha=0.4,linestyle='--')
ax.set_xticks(np.arange(1,32,2))
ax.set_xticklabels([str(i) for i in range(1, 32,2)],fontsize=16)
ax.set_yticks(np.arange(0,31,5))
ax.set_yticklabels([str(i) for i in range(0, 31,5)],fontsize=16)



for t in unique_times:
    idx = np.where( time_dates == t )
    N=np.ones(len(idx[0]))*n

    lower.append(np.where( val.data[idx[0]]*0.3 <0.02 )[0].shape[0])
    
    # plt.plot(np.ones(len(idx[0])),val[idx],'r')
    plt.plot(N,err[idx],'b',marker='o',linewidth=0,markersize=3)
    plt.plot(N,val[idx],'r',marker='o',linewidth=0,markersize=3)
    n+=1
    a12.append(np.sum(val[idx]>=11))
    data_len.append(len(idx[0]))
ax.legend(['obs_error','obs_values'])













    # idx = np.where( (time_dates >= dt.datetime(2025,5,1)) & (time_dates <= dt.datetime(2025,5,8) ))


#     grid_val = griddata(
#     (lon[idx], lat[idx]), val[idx],
#     (grid_x, grid_y), method='linear'  # or 'cubic', 'nearest'
#     )

#     fig = plt.figure(figsize=(10, 6))
#     ax = plt.axes(projection=ccrs.PlateCarree())
#     ax.set_title(f'Phyt absolute value- (Less than 0.02^2)', fontsize=14)
#     # ax.set_title(f'Phyt Observation Error', fontsize=14)

#     ax.set_extent([105, 165, 15, 50], crs=ccrs.PlateCarree())

#     # 해안선 등 추가
#     ax.coastlines(resolution='10m', linewidth=1)
#     ax.add_feature(cfeature.BORDERS, linestyle=':')
#     ax.add_feature(cfeature.LAND, facecolor='lightgray')
#     ax.add_feature(cfeature.OCEAN, facecolor='#ffffff')
#     ax.gridlines(draw_labels=True, linestyle='--')

#     # 관측 산점도
#     val_min=np.nanmin(val[idx])
#     idV=np.nanargmin(val[idx])
    
#     lonV=lon[idV]
#     latV=lat[idV]

#     ax.text(
#     0.01, 0.99,  # (x, y): 좌측 상단 기준 (axes fraction)
#     f'Min: {val_min:.2f} mmol/m^3\n'+f'Min lon: {lonV:.2f}\n'+f'Min lat: {latV:.2f}' ,
#     transform=ax.transAxes,  # 축 좌표계 기준
#     fontsize=12,
#     fontweight='bold',
#     va='top',
#     ha='left',
#     bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')  # 배경 박스 옵션
# )
#     sc = ax.scatter(lon[idx], lat[idx], c=val[idx], cmap='jet', s=10,vmin=.01, vmax=.02,
#                     edgecolor=None, transform=ccrs.PlateCarree())

#     # 컬러바
#     cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.06, shrink=.8)
#     cbar.set_label('(mmol m-3)', fontsize=12)

#     plt.tight_layout()
#     # plt.savefig(f'D:/shjo/ROMS_OUTS/E_test/Case_03E/figs/obs_error/phyto_err_{t.strftime("%Y%m%d")}', dpi=200)

#     plt.show()
    
    
    
    
    
    
    