

# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import create_I as cn
import utils as tl
from io_utils import collect_time_info
print("[💣] STARTING FULL EXECUTION...")
# --- [01] Load configuration and input metadata ---
cfg  = tl.parse_config("./config.yaml")
grd  = tl.load_roms_grid(cfg.grdname)
ogcm = tl.load_ogcm_metadata(cfg.ogcm_name, cfg.ogcm_var_name)


# --- [02] Time index matching and relative time calculation ---
tinfo = collect_time_info(cfg.ogcm_name, cfg.ogcm_var_name['time'], str(cfg.initdate))
_, idt, _ = tinfo[0]
relative_time = tl.compute_relative_time(ogcm.time[idt], ogcm.time_unit, cfg.time_ref)

# --- [03] Create initial NetCDF file ---
print(f"--- [03] Creating initial NetCDF file: {cfg.ininame} ---")
status = cn.create_ini(cfg, grd, relative_time, ncFormat=cfg.ncformat, bio_model="Fennel")
if status:
    print(f"--- [!ERROR] Failed to creating file {cfg.ininame} ---")
    raise
print(f"--- Created file: {cfg.ininame} ---")

# --- [04] Crop OGCM domain and prepare remap weights ---
lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)

if cfg.calc_weight:
    print(f"--- [04] Calculating weight: {cfg.weight_file} ---")
    status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
    
    if status:
        print(f"--- [!ERROR] Failed to generate remap weights: {cfg.weight_file} ---")
        raise  
else: 
    print(f"--- [04] Use existing wght file {cfg.weight_file} ---")


# --- [05] Load OGCM raw fields (zeta, temp, salt, u, v) ---
print("--- [05] Loading OGCM data ---")
with Dataset(cfg.ogcm_name, maskandscale=False) as nc_raw:
    nc = tl.MaskedNetCDF(nc_raw)

    zeta = nc.get(cfg.ogcm_var_name['zeta'], idt, idy, idx)
    temp = nc.get(cfg.ogcm_var_name['temperature'], idt, slice(None), idy, idx)
    salt = nc.get(cfg.ogcm_var_name['salinity'], idt, slice(None), idy, idx)
    u    = nc.get(cfg.ogcm_var_name['u'], idt, slice(None), idy, idx)
    v    = nc.get(cfg.ogcm_var_name['v'], idt, slice(None), idy, idx)
    ubar = tl.depth_average(u, ogcm.depth)
    vbar = tl.depth_average(v, ogcm.depth)

field = tl.ConfigObject(
    zeta = zeta,
    ubar = ubar,
    vbar = vbar,
    temp = temp,
    salt = salt,
    u    = u,
    v    = v
)


# --- [06] Load and apply remap weights to all fields ---
print("--- [06] Remapping ---")
with Dataset(cfg.weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]

for varname in vars(field) :
    var_src = getattr(field, varname)
    remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
    #remapped = np.nan_to_num(remapped, nan=0.0)
    setattr(field, varname, remapped)

# Debugging
print("[dbg] Check coastal NaN values")
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")

dbg_val=field.temp[0,...]

filled_dbg = tl.fast_flood_2d(dbg_val)
filled_dbg = tl.fast_flood_3d(field.temp)

tmp=np.ones_like(dbg_val)
tmp[grd.mask==0]=np.nan
tmp[dbg_val==dbg_val]=np.nan

plt.figure()
plt.pcolor(grd.lon,grd.lat,field.temp[-2])
plt.clim([1.2,1.7])
plt.colorbar()
plt.savefig("/mnt/c/Users/shjo9/Bridge/temp_bottom.png")
plt.close

plt.figure()
plt.pcolor(grd.lon,grd.lat,filled_dbg[-2])
plt.clim([1.2,1.7])
plt.colorbar()
plt.savefig("/mnt/c/Users/shjo9/Bridge/temp__filled_bottom.png")
plt.close



raise






# --- [07] Rotate vectors to model grid and convert to u/v points ---
print("--- [07] Rotate vectors ---")
u_rot, v_rot       = tl.rotate_vector_euler(field.u,    field.v,    grd.angle, to_geo=False)
ubar_rot, vbar_rot = tl.rotate_vector_euler(field.ubar, field.vbar, grd.angle, to_geo=False)

setattr(field, 'u',    tl.rho2uv(u_rot,'u'))
setattr(field, 'v',    tl.rho2uv(v_rot,'v'))
setattr(field, 'ubar', tl.rho2uv(ubar_rot,'u'))
setattr(field, 'vbar', tl.rho2uv(vbar_rot,'v'))

# --- [08] Vertical interpolation from z-level to sigma-level ---
print("--- [08] Vertical interpolation from z-level to sigma-level ---")
zlevs_args = (
    cfg.vertical.vtransform,
    cfg.vertical.vstretching,
    cfg.vertical.theta_s,
    cfg.vertical.theta_b,
    cfg.vertical.tcline,
    cfg.vertical.layer_n,
)

zr = tl.zlevs(*zlevs_args, 1, grd.topo, field.zeta)
zu, zv = tl.rho2uv(zr,'u'), tl.rho2uv(zr,'v')

Z=np.zeros(len(ogcm.depth)+2)
Z[0]=100;Z[1:-1]=-ogcm.depth;Z[-1]=-100000

for var, zgrid in zip(['temp', 'salt', 'u', 'v'], [zr, zr, zu, zv]):
    val = getattr(field, var)
    padded = np.vstack((val[0:1], val, val[-1:]))
    flipped = np.flip(padded, axis=0)
    remapped = tl.ztosigma(flipped, zgrid, np.flipud(Z))
    setattr(field, var, remapped)


# --- [09] Volume conservation and barotropic velocity correction ---
print(" --- [09] Volume conservation and barotropic velocity correction ---")
zw = tl.zlevs(*zlevs_args, 5, grd.topo, field.zeta)
dzr = zw[1:, :, :] - zw[:-1, :, :]
dzu, dzv = tl.rho2uv(dzr,'u'), tl.rho2uv(dzr,'v')

u_consv, v_consv, ubar_new, vbar_new = \
    tl.conserve_and_recompute_barotropic(
        field.u, field.v,
        field.ubar, field.vbar,
        dzu, dzv
    )

setattr(field, 'u', u_consv)
setattr(field, 'v', v_consv)
setattr(field, 'ubar', ubar_new)
setattr(field, 'vbar', vbar_new)

# --- [10] Write all remapped variables to ini.nc ---
print(f"--- [10] Write all remapped variables to {cfg.ininame} ---")
with Dataset(cfg.ininame, mode='a') as nc:
    for var in ['zeta', 'temp', 'salt', 'u', 'v', 'ubar', 'vbar']:
        nc[var][0] = getattr(field, var)


print("--- [DONE] ---")














