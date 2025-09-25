from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

pth = "/data/shjo/data/nifs01"

ncG_src=Dataset(pth+'/NWP4_grd_314_10m.nc')
ncG_dst=Dataset(pth+'/source/NWP12_grd_NWP4.nc','a')

lon_src,lat_src=ncG_src["lon_rho"][:],ncG_src["lat_rho"][:]
lon_dst,lat_dst=ncG_dst["lon_rho"][:],ncG_dst["lat_rho"][:]

topo_src=ncG_src["h"][:].data
#topo_dst=ncG_dst["h"][:].data

#print(topo_src.shape)
#print(topo_dst.shape)

topo_intrp = griddata( (lon_src.flatten(), lat_src.flatten()), topo_src.flatten(), (lon_dst,lat_dst) )

#print(topo_intrp.shape)
ncG_dst["h"][:]=topo_intrp

ncG_dst.close()
ncG_src.close()













