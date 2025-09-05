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

pthG_src="/home/shjo/data/roms_inputs/grd/lte/NWP4_grd_3_10m_LP.nc"
pthG_dst="/home/shjo/data/roms_inputs/grd/lte/NWP12_grd_NWP4.nc"
pthI_src="/home/shjo/data/lab/NWP4_ini_DA_8401.nc"
pthI_dst="/home/shjo/data/lab/NWP4_ini_DA_8401.nc"

weight_file="./test_wght.nc"

# Read grid
ncG_src=Dataset(pthG_src)
ncG_dst=Dataset(pthG_dst)
ncI_src=Dataset(pthI_src)
ncI_dst=Dataset(pthI_dst)

lon_src,lat_src=ncG_src['lon_rho'][:], ncG_src['lat_rho'][:]
lon_dst,lat_dst=ncG_dst['lon_rho'][:], ncG_dst['lat_rho'][:]

topo_src = ncG_src["h"][:]
topo_dst = ncG_dst["h"][:]
mask_dst = ncG_dst["mask_rho"][:]

angle_src = ncG_src["angle"][:]
angle_dst = ncG_dst["angle"][:]

vtrs=float(ncI_src["Vtransform"][:].data)
vstr=float(ncI_src["Vstretching"][:].data)
theta_s=float(ncI_src["theta_s"][:].data)
theta_b=float(ncI_src["theta_b"][:].data)
theta_b=0.1
tcline=float(ncI_src["Tcline"][:].data)
nlayer=len(ncI_src["sc_r"][:])

zr_src=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer,1,topo_src,np.zeros_like(topo_src))
zr_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,36,1,topo_dst,np.zeros_like(topo_dst))


# --- crop domain ---
lon_crop, lat_crop, idx, idy = ut.crop_to_model_domain(lat_dst, lon_dst, lat_src, lon_src)

# --- calc weight ---
status = ut.build_bilinear_regridder(lon_src, lat_src, lon_dst, lat_dst, weight_file, reuse=False)
if status:
    raise RuntimeError(f"Failed to generate remap weights: {weight_file}")

# --- Load & calc increments ---
factor=1
zeta_diff=(ncI_src["zeta"][0]-ncI_src["zeta"][1])*factor
temp_diff=(ncI_src["temp"][0]-ncI_src["temp"][1])*factor
salt_diff=(ncI_src["salt"][0]-ncI_src["salt"][1])*factor
u_diff=(ncI_src["u"][0]-ncI_src["u"][1])*factor
v_diff=(ncI_src["v"][0]-ncI_src["v"][1])*factor
ubar_diff=(ncI_src["ubar"][0]-ncI_src["ubar"][1])*factor
vbar_diff=(ncI_src["vbar"][0]-ncI_src["vbar"][1])*factor

u_diff=pu.uv2rho_rutgers_safenan(u_diff,"u")
v_diff=pu.uv2rho_rutgers_safenan(v_diff,"v")
ubar_diff=pu.uv2rho_rutgers_safenan(ubar_diff,"u")
vbar_diff=pu.uv2rho_rutgers_safenan(vbar_diff,"v")

# --- define field ---
field = ut.ConfigObject(zeta=zeta_diff, temp=temp_diff, salt=salt_diff,\
        u=u_diff,v=v_diff,ubar=ubar_diff,vbar=vbar_diff,angle=angle_src)

# --- Load and apply remap weights to all fields ---
with Dataset(weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]
for varname in vars(field):
    var_src = getattr(field, varname)
    remapped = ut.remap_variable(var_src, row, col, S, lon_dst.shape, method="coo")
    setattr(field, varname, remapped)

# --- Horizontal flood (all fields) ---
for var in vars(field):
    val = getattr(field, var)
    val_flooded = ut.flood_horizontal(val, lon_dst, lat_dst, method="griddata")
    setattr(field, var, val_flooded)
#    done("flood_h", time.time()-t0)

# --- interp zr ---
zr_src_remapped = ut.remap_variable(zr_src, row, col, S, lon_dst.shape, method="coo")
zr_src_flooded = ut.flood_horizontal(zr_src_remapped, lon_dst, lat_dst, method="griddata")

# --- Mask land to 0 ---
for varname in vars(field):
    var = getattr(field, varname)
    if var.ndim==2:
        var[mask_dst == 0] = 0.0
    else:
        var[:, mask_dst == 0] = 0.0
    setattr(field, varname, var)

# --- vertical interpolation ---
for varname in vars(field):
    var = getattr(field,varname)
    if var.ndim==2:
        continue
    print(varname)
    var = getattr(field, varname)
    var_zinterp = pu.vertical_interp_to_ZR(zr_src_flooded, var, zr_dst,
                                         n_jobs=-1,
                                         dedup="mean",
                                         extrap_mode="padding")
    setattr(field, varname, var_zinterp)

# --- rotate vector ---
# 1. parent grid -> geo (rho-point 회전)
u,v=getattr(field,"u"), getattr(field,"v")
angle_re=getattr(field,"angle")
Urho, Vrho = ut.rotate_vector_euler(u, v, angle_re, to_geo=True)

# 2. geo -> target grid (rho-point 회전)
Ugeo, Vgeo = ut.rotate_vector_euler(Urho, Vrho, angle_dst, to_geo=False)

# 3. staggered 분배
setattr(field,"u",ut.rho2uv(Ugeo, pos="u"))
setattr(field,"v",ut.rho2uv(Vgeo, pos="v"))




