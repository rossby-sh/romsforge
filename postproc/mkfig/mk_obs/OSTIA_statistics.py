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

Spth='D:/shjo/LTE/OSTIA_SST_2023-2024.nc'
ncX=xr.open_dataset(Spth).analysed_sst.loc[dict(time=slice('2023-04-01','2023-04-30'))] - 273.15
SST=ncX.mean(dim='time').values
lonO_,latO_=ncX.longitude.values,ncX.latitude.values
lonO,latO=np.meshgrid(lonO_,latO_)

Gpth='D:/shjo/ROMS_inputs/NWP4_grd_3_10m_LP.nc'
ncG=Dataset(Gpth)
lon_rho,lat_rho=ncG['lon_rho'][:],ncG['lat_rho'][:]

zi = griddata((lonO.flatten(),latO.flatten()), SST.flatten(), (lon_rho,lat_rho), method='linear')

nrpth='C:/Users/ust21/shjo/LTE/LTE_OL_para_apr.nc'
# nrpth='C:/Users/ust21/shjo/LTE/LTE_control2_apr.nc'

nc=Dataset(nrpth)
SSTM=nc['temp'][0,-1]



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
fac=lnum_ECS


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
prmse = (1 - (0.79445-rmse)/rmse*100)
print('ES %rmse: '+str(prmse)[:3])
prmse = (1 - ( 0.56072-rmse)/rmse*100)
print('YS %rmse: '+str(prmse)[:3])
prmse = (1 - (0.79445-rmse)/rmse*100)
print('ECS %rmse: '+str(prmse)[:3])



