from netCDF4 import Dataset, num2date, date2num
import datetime as dt

pth="hycom_korea_20260311.nc"

nc=Dataset(pth,"a")

time=nc["time"][:]

time_re=time+47

print(num2date(time_re,"days since 2000-1-1"))

nc["time"][:]=time_re

nc.close()






