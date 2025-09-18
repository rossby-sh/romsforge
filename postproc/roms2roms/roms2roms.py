# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
#import create_I as cn
import utils as ut
import post_utils as pu

pthG_pre="/home/shjo/data/roms_inputs/grd/lte/NWP4_grd_3_10m_LP.nc"
pthG_pst="/home/shjo/data/roms_inputs/grd/lte/NWP12_grd_NWP4.nc"
pthI="/home/shjo/data/lab/NWP4_ini_DA_8401.nc"

weight_file="./test_wght.nc"
#vdepth = -np.array([0.5, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000,\
#          1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 5500])

ncG_pre=Dataset(pthG_pre)
ncG_pst=Dataset(pthG_pst)
ncI=Dataset(pthI)

lon_pre,lat_pre=ncG_pre['lon_rho'][:], ncG_pre['lat_rho'][:]
lon_pst,lat_pst=ncG_pst['lon_rho'][:], ncG_pst['lat_rho'][:]

topo_pre = ncG_pre["h"][:]
topo_pst = ncG_pst["h"][:]
mask_pst = ncG_pst["mask_rho"][:]

vtrs=float(ncI["Vtransform"][:].data)
vstr=float(ncI["Vstretching"][:].data)
theta_s=float(ncI["theta_s"][:].data)
theta_b=float(ncI["theta_b"][:].data)
theta_b=0.1
tcline=float(ncI["Tcline"][:].data)
nlayer=len(ncI["sc_r"][:])

zr_pre=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer,1,topo_pre,np.zeros_like(topo_pre))
zr_pst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,36,1,topo_pst,np.zeros_like(topo_pst))


# crop domain
lon_crop, lat_crop, idx, idy = ut.crop_to_model_domain(lat_pst, lon_pst, lat_pre, lon_pre)
status = ut.build_bilinear_regridder(lon_pre, lat_pre, lon_pst, lat_pst, weight_file, reuse=False)
if status:
    raise RuntimeError(f"Failed to generate remap weights: {weight_file}")

# --- Load parent initial data ---
idt=0
with Dataset(pthI, maskandscale=True) as nc_raw:
    nc = ut.MaskedNetCDF(nc_raw)
    zeta = nc.get('zeta',        idt, slice(None))
    temp = nc.get('temp',        idt, slice(None))
    salt = nc.get('salt',        idt, slice(None))
    u    = nc.get('u',           idt, slice(None))
    v    = nc.get('v',           idt, slice(None))
    ubar = nc.get("ubar",        idt, slice(None))
    vbar = nc.get("vbar",        idt, slice(None))
field = ut.ConfigObject(zeta=zeta, temp=temp, salt=salt)

print(getattr(field,"temp").shape)

# [06] Load and apply remap weights to all fields
with Dataset(weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]
for varname in vars(field):
    var_src = getattr(field, varname)
    remapped = ut.remap_variable(var_src, row, col, S, lon_pst.shape, method="coo")
    setattr(field, varname, remapped)

# [07] Horizontal flood (all fields)
for var in ['temp', 'salt','zeta']:
    val = getattr(field, var)
    val_flooded = ut.flood_horizontal(val, lon_pre, lat_pst, method="griddata")
    setattr(field, var, val_flooded)
#    done("flood_h", time.time()-t0)

# [09] Mask land to 0
for varname in ['zeta']:
    var = getattr(field, varname)
    var[mask_pst == 0] = 0.0
    setattr(field, varname, var)
for varname in ['temp', 'salt']:
    var = getattr(field, varname)
    var[:, mask_pst == 0] = 0.0
    setattr(field, varname, var)

zr_pre_remapped = ut.remap_variable(zr_pre, row, col, S, lon_pst.shape, method="coo")
zr_pre_flooded = ut.flood_horizontal(zr_pre_remapped, lon_pre, lat_pst, method="griddata")

for varname in ["temp", "salt"]:
    print(varname)
    var = getattr(field, varname)
    var_zinterp = pu.vertical_interp_to_ZR(zr_pre_flooded, var, zr_pst,
                                         n_jobs=-1,
                                         dedup="mean",
                                         extrap_mode="leading")
    setattr(field, varname, var_zinterp)







