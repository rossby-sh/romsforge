# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 10:56:45 2025

@author: ust21
"""

import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import numpy as np
import sys
sys.path.append('C:/Users/ust21/shjo/others/will_delet/myROMS/prc_src/utils/')
from netCDF4 import Dataset, num2date
import ROMS_utils02 as ru
import cmocean
import xarray as xr
from scipy.interpolate import griddata
from copy import deepcopy
import datetime
import xarray as xr

# 파일 열기
pth = 'D:/shjo/tmp/2023/'
flist=[pth+i for i in os.listdir(pth) if i.endswith('.nc')]

ref_date = np.datetime64('2000-01-01')
datasets=[]

for f in flist:
    ds = xr.open_dataset(f)
    ds = ds.expand_dims('time')

    # 날짜 추출
    fname = os.path.basename(f)
    date_str = fname.split('.')[1]  # 예: '20241231'
    dt = datetime.datetime.strptime(date_str, "%Y%m%d")
    np_time = (np.datetime64(dt) - ref_date) / np.timedelta64(1, 'D')  # float days

    ds['time'] = ('time', [np_time])
    ds['time'].attrs['units'] = 'days since 2000-01-01 00:00:00'
    datasets.append(ds)


ds_combined = xr.concat(datasets, dim='time').sel(
    lon=slice(100,170), lat=slice(55,5))
val=ds_combined.chlor_a.mean(dim='time',skipna=True).values#/(0.02*6.625*12)

lon, lat = ds_combined.lon.values, ds_combined.lat.values

lonO,latO=np.meshgrid(lon,lat)




Gpth='D:/shjo/ROMS_inputs/NWP4_grd_3_10m_LP.nc'
ncG=Dataset(Gpth)
lon_rho,lat_rho=ncG['lon_rho'][:],ncG['lat_rho'][:]

zi = griddata((lonO.flatten(),latO.flatten()), val.flatten(), (lon_rho,lat_rho), method='linear')

nrpth='C:/Users/ust21/shjo/LTE/LTE_OL_para_apr.nc'
# nrpth='C:/Users/ust21/shjo/LTE/LTE_control2_apr.nc'

nc=Dataset(nrpth)
SSTM=nc['chlorophyll'][0,-1]



# Regionnal 

# lnum_ES = np.where((lon_rho > 127.0) & (lon_rho < 143.0) &
#                    (lat_rho > 33.0 ) & (lat_rho < 49.0 ))
# lnum_YS = np.where((lon_rho > 117.0) & (lon_rho < 127.0) &
#                    (lat_rho > 33.0 ) & (lat_rho < 42.0 ))
# lnum_ECS = np.where((lon_rho > 119.0) & (lon_rho < 126.0) &
#                     (lat_rho > 26.0 ) & (lat_rho < 33.0 ))
# lnum_KOE = np.where((lon_rho > 141.0) & (lon_rho < 155.0) &
#                     (lat_rho > 33.0 ) & (lat_rho < 43.0 ))


lnum_ES = np.where((lon_rho <= 127.0) | (lon_rho >= 143.0) |
                   (lat_rho <= 33.0 ) | (lat_rho >= 49.0))

lnum_YS = np.where((lon_rho <= 117.0) | (lon_rho >= 127.0) |
                   (lat_rho <= 33.0 ) | (lat_rho >= 42.0))

lnum_ECS = np.where((lon_rho <= 119.0) | (lon_rho >= 126.0) |
                    (lat_rho <= 26.0 ) | (lat_rho >= 33.0))

lnum_KOE = np.where((lon_rho <= 141.0) | (lon_rho >= 155.0) |
                    (lat_rho <= 33.0 ) | (lat_rho >= 43.0))


Model= deepcopy(SSTM)
obs  = deepcopy(zi)
fac=lnum_ES

Model[fac]=np.nan
obs[fac]=np.nan
obs[Model!=Model]=np.nan
Model[obs!=obs]=np.nan


# calc bias
bias = Model - obs
mean_bias = np.nanmean(bias)
print('mean_bias: '+str(mean_bias))

# calc nRMSE
rmse = np.sqrt(np.nanmean(bias**2))
nrmse1 = rmse / np.nanstd(obs.data)  # % 단위
print('nRMSE01: '+str(nrmse1))

# obs_range = np.nanmax(obs) - np.nanmin(obs)
# nrmse2 = rmse / obs_range * 100  # % 단위
# print('nRMSE02: '+str(nrmse2))

# calc Corr
corr_matrix = np.corrcoef(Model[Model==Model], obs[Model==Model])
corr = np.nanmean(corr_matrix)  # 상관계수
print('Corr: '+str(corr)[:6])

# calc nRMSE %
prmse = (1 - (2.864612-rmse)/rmse*100)
print('ES %rmse: '+str(prmse)[:3])
prmse = (1 - (1.9772498671698384-rmse)/rmse*100)
print('YS %rmse: '+str(prmse)[:3])
prmse = (1 - (0.871021-rmse)/rmse*100)
print('ECS %rmse: '+str(prmse)[:3])



