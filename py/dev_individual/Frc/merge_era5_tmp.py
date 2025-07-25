### 
#
# Auther: Seonghyun Jo
# Description: This code merges downloaded ERA5 data from C3S and slice it to NWP domain
#
###

import xarray as xr

rpth='/home/shjo/DATA/era5_tmp2/'
wnpth='/data/share/DATA/PROC/ERA5/ERA5_merged_20250531-20250701.nc'

Accum = xr.open_mfdataset(rpth+'data_stream-oper_stepType-accum*').loc[dict(latitude=slice(60,5),longitude=slice(100,170))] 
Instant = xr.open_mfdataset(rpth+'data_stream-oper_stepType-instant*').loc[dict(latitude=slice(60,5),longitude=slice(100,170))] 

print('=== Reading Data ===')
ERA5 = xr.merge([Accum,Instant],compat='override')

#print('=== ERA5 valid_time ===')
#print(ERA5.valid_time.values)

print('=== Writing to '+wnpth+' ===')

ERA5.to_netcdf(wnpth)

