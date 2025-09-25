from netCDF4 import Dataset

pth="/data/shjo/data/nifs01/nifs01_8km_river.nc"

nc=Dataset(pth,"a")
nc["river_time"].cycle_length = 365.25
nc.close()



