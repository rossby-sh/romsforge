# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
from netCDF4 import Dataset, num2date, date2num
import time
import datetime as dt

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_I as cn
import utils as tl
from io_utils import collect_time_info_legacy

# logging
from log_utils2 import configure, step, capture_warnings, info, note, plus, bar, warn_line, done, ellipsis
configure(width=80, show_sections=False, color_mode='auto')

# --- main --------------------------------------------------------------------
start0 = time.time()
bar("Initial Condition Build (ini)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg  = tl.parse_config("./config_clm.yaml")
    grd  = tl.load_roms_grid(cfg.grdname)
    ogcm = tl.load_ogcm_metadata(cfg.bio_ogcm_name, cfg.bio_ogcm_var_name)
    info(f"grid={cfg.grdname}")
    info(f"ini_out={cfg.ininame}")
    info(f"wght_out={cfg.bio_weight_file}")
    info(f"ogcm={cfg.bio_ogcm_name}")

with step("[02] Load climatology file"):
    clim_file = cfg.bio_ogcm_name
    assert os.path.exists(clim_file), f"Climatology file not found: {clim_file}"
    relative_time = date2num(cfg.initdate,cfg.time_ref)
    print(relative_time)
    info(f"Loaded climatology file: {os.path.basename(clim_file)}")

# [03] Create initial NetCDF file
#with step("[03] Create initial NetCDF"):
#    status = cn.create_ini__(cfg, grd, relative_time, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
#    if status:
#        raise RuntimeError(f"Failed creating file {cfg.ininame}")

# [04] Crop OGCM domain and prepare remap weights
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.bio_weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.bio_weight_file}")
    else:
        info(f"Use existing wght file {cfg.bio_weight_file}")

nc_clim = tl.MaskedNetCDF(Dataset(cfg.bio_ogcm_name, mode="r", maskandscale=True))
with step("[05] Load OGCM fields", ts="Load OGCM"):
    NO3  = nc_clim.get(cfg.bio_ogcm_var_name['NO3'], 0, slice(None), idy, idx)
    phyt = nc_clim.get(cfg.bio_ogcm_var_name['phyt'], 0, slice(None), idy, idx)
    chlo = nc_clim.get(cfg.bio_ogcm_var_name['chl'], 0, slice(None), idy, idx)
    PO4 = nc_clim.get(cfg.bio_ogcm_var_name['PO4'], 0, slice(None), idy, idx)
    oxygen = nc_clim.get(cfg.bio_ogcm_var_name['oxygen'], 0, slice(None), idy, idx)
    zoop = phyt*0.3

    alkalinity = np.ones_like(zoop) * 2350
    TIC = np.ones_like(zoop) * 2100
    NH4 = np.ones_like(zoop) * 0.01

    SdeC = np.ones_like(zoop) * 0.04
    LdeC = np.ones_like(zoop) * 0.04
    RdeC = np.ones_like(zoop) * 0.04
    SdeN = np.ones_like(zoop) * 0.04
    RdeN = np.ones_like(zoop) * 0.04

    field = tl.ConfigObject(NO3=NO3, phytoplankton=phyt, chlorophyll=chlo, PO4=PO4, oxygen=oxygen, zooplankton=zoop,alkalinity=alkalinity, TIC=TIC, NH4=NH4, SdetritusC=SdeC, LdetritusC=LdeC,RdetritusC=RdeC, SdetritusN=SdeN, RdetritusN=RdeN)

# [06] Load and apply remap weights to all fields
with step("[06] Remap (weights)"):
    with Dataset(cfg.bio_weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S   = nc.variables["S"][:]
    t0 = time.time()
    for varname in vars(field):
        var_src = getattr(field, varname)
        with capture_warnings(tag="remap"):
            remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
        setattr(field, varname, remapped)
#    done("remap", time.time()-t0)

note(f"Flood method: {cfg.flood_method_for_ini}")
# [07] Horizontal flood (all fields)
with step("[07] Flood: horizontal"):
    t0 = time.time()
    for var in vars(field):
        val = getattr(field, var)
        val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_ini)
        setattr(field, var, val_flooded)
#    done("flood_h", time.time()-t0)

# [08] Vertical flood
with step("[08] Flood: vertical"):
    t0 = time.time()
    for var in vars(field) :
        val = getattr(field, var)
        with capture_warnings(tag="vflood"):
            val_flooded = tl.flood_vertical_vectorized(val, grd.mask_rho, spval=-1e10)
            # Alternative:
            # val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask), spval=-1e10)
        setattr(field, var, val_flooded)
#    done("flood_v", time.time()-t0)

# [09] Mask land to 0
with step("[10] Mask & rotate", ts=f"None"):
    for var in vars(field):
        arr = getattr(field, var)
        if arr.ndim == 2:
            arr[grd.mask_rho == 0] = 0.0
        else:
            arr[..., grd.mask_rho == 0] = 0.0
        setattr(field, var, arr)

# [11] Vertical interpolation from z-level to sigma-level
with step("[11] z→σ interpolation"):
    zlevs_args = (
        cfg.vertical.vtransform,
        cfg.vertical.vstretching,
        cfg.vertical.theta_s,
        cfg.vertical.theta_b,
        cfg.vertical.tcline,
        cfg.vertical.layer_n,
    )
    zr = tl.zlevs(*zlevs_args, 1, grd.topo, np.zeros_like(grd.topo))

    Z = np.zeros(len(ogcm.depth)+2)
    Z[0] = 100; Z[1:-1] = -np.abs(ogcm.depth); Z[-1] = -100000
    Zf = np.flipud(Z)
    for var, zgrid in zip(["NO3","PO4","NH4","TIC","oxygen","chlorophyll","phytoplankton","zooplankton","alkalinity","SdetritusC","LdetritusC","RdetritusC","SdetritusN","RdetritusN"], [zr, zr, zr, zr,zr,zr,zr,zr,zr,zr,zr,zr,zr,zr]):
        val = getattr(field, var)
        padded = np.vstack((val[0:1], val, val[-1:]))
        flipped = np.flip(padded, axis=0)
        with capture_warnings(tag="z2sigma"):
            remapped = tl.ztosigma_numba(flipped, zgrid, Zf)
        setattr(field, var, remapped)

bar("Summary")
print(f"Total elapsed: {time.time() - start0:.3f}s")
# [12] Write all remapped variables to bry.nc
with step("[12] Write variables", out=ellipsis(cfg.ininame)):
    with Dataset(cfg.ininame, 'a') as nc:
        for varname in vars(field):
            nc[varname][0] = getattr(field,varname)

bar("Summary")
