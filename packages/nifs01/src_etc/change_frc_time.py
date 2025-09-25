
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import datetime as dt
pth="/data/shjo/data/nifs01/source/nifs_frc_inst_19900101-20250731.nc"


nc=Dataset(pth,'a')

time=nc['wind_time'][:].data
time_unit=nc["wind_time"].units
print(time[:5])

diff=date2num(dt.datetime(2000,1,1),"seconds since 2000-1-1") - date2num(dt.datetime(2000,1,1),time_unit)
print(diff)
print(num2date(time[:5]+diff,"seconds since 2000-1-1"))
print(time[:5]+diff)

T=time+diff

nc["wind_time"][:]=T
nc["tair_time"][:]=T
nc["qair_time"][:]=T
nc["cloud_time"][:]=T
nc["sst_time"][:]=T
nc["dqdsst_time"][:]=T
nc["srf_time"][:]=T
nc["lrf_time"][:]=T
nc["rain_time"][:]=T
nc["pair_time"][:]=T

nc["wind_time"].units="seconds since 2000-1-1"
nc["tair_time"].units="seconds since 2000-1-1"
nc["qair_time"].units="seconds since 2000-1-1"
nc["cloud_time"].units="seconds since 2000-1-1"
nc["sst_time"].units="seconds since 2000-1-1"
nc["dqdsst_time"].units="seconds since 2000-1-1"
nc["srf_time"].units="seconds since 2000-1-1"
nc["lrf_time"].units="seconds since 2000-1-1"
nc["rain_time"].units="seconds since 2000-1-1"
nc["pair_time"].units="seconds since 2000-1-1"

nc.close()

