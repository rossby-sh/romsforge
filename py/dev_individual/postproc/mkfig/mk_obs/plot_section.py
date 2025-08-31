# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 19:02:39 2025

@author: ust21
"""


import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import numpy as np

from netCDF4 import Dataset, num2date
# import cmocean
import sys
sys.path.append('C:/Users/ust21/shjo/others/will_delet/myROMS/prc_src/utils')
import ROMS_utils01 as ru


Gpth='D:/shjo/ROMS_inputs/NWP12_grd_NWP4.nc'
rpth='D:/shjo/LTE/ECCO_MIRES/avg/NWP12_avg_8449_0008.nc'

ncG=Dataset(Gpth)
ncA=Dataset(rpth)


vname = 'temp'
lon_rng=[156]
lat_rng=[15,20.5]

[X,Z,VAR] = ru.get_section(ncG,ncA,vname,lon_rng,lat_rng,tindx=0)


plt.pcolor(X,Z,VAR,cmap=plt.get_cmap('Spectral_r',27))
plt.colorbar()
plt.ylim([-1000,0])


vname = 'temp'
lon_rng=[150,157]
lat_rng=[15.2]

[X,Z,VAR] = ru.get_section(ncG,ncA,vname,lon_rng,lat_rng,tindx=0)


plt.pcolor(X,Z,VAR,cmap=plt.get_cmap('Spectral_r',27))
plt.colorbar()
plt.ylim([-1000,0])


import xarray as xr

A=xr.open_dataset(rpth).temp[-1,-1].loc[dict(eta_rho=slice(0,20),xi_rho=slice(510,600))]

A.plot(cmap=plt.get_cmap('Spectral_r',27))

A.min()





















