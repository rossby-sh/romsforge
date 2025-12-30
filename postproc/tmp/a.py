from netCDF4 import Dataset



nc=Dataset("/home/shjo/github/romsforge/postproc/tmp/nifs5km_avg_9405_0001.nc","a")




nc["ocean_time"][0] = 9404*86400

nc.close()



