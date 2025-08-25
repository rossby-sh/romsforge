### 
#
# Auther: Seonghyun Jo
# Description: This code merges downloaded ERA5 data from C3S and slice it to NWP domain
#
###

import xarray as xr

#rpth='/data/share/DATA/RAW/ERA5/'
rpth='/data/share/DATA/PROC/ERA5/monthly/'
#wnpth='/data/share/DATA/PROC/ERA5/ERA5_221231-230228.nc'
wnpth='/data/share/DATA/PROC/ERA5/ERA5_monthly_30Y.nc'

Accum = xr.open_mfdataset('/data/share/DATA/PROC/ERA5/monthly/data_stream-moda_stepType-avgad.nc').loc[dict(latitude=slice(60,5),longitude=slice(100,175))] 
Instant = xr.open_mfdataset('/data/share/DATA/PROC/ERA5/monthly/data_stream-moda_stepType-avgua.nc').loc[dict(latitude=slice(60,5),longitude=slice(100,175))] 

print('=== Reading Data ===')
ERA5 = xr.merge([Accum,Instant],compat='override')

#print('=== ERA5 valid_time ===')
#print(ERA5.valid_time.values)

print('=== Writing to '+wnpth+' ===')

ERA5.to_netcdf(wnpth)

