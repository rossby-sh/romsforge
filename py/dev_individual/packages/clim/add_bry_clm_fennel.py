# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
from collections import defaultdict
import time

# Append libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_B as cn
import utils as tl
from io_utils import collect_time_info
from log_utils2 import configure, step, capture_warnings, info, note, plus, bar, warn_line, done, ellipsis
configure(width=80, show_sections=False)  # ruler 끄기

# === pretty log helpers =======================================================
# --- main --------------------------------------------------------------------
start1 = time.time()
bar("Boundary Build")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg = tl.parse_config("./config_single.yaml")
    grd = tl.load_roms_grid(cfg.grdname)
    ogcm = tl.load_ogcm_metadata(cfg.bio_ogcm_name, cfg.bio_ogcm_var_name)
    info(f"src=OGCM grid={os.path.basename(cfg.grdname)}")
    info(f"out={ellipsis(cfg.bryname)} wght={ellipsis(cfg.weight_file)}")

# --- [02] Load climatology file and collect time info ---
with step("[02] Load climatology file"):
    clim_file = cfg.bio_ogcm_name
    assert os.path.exists(clim_file), f"Climatology file not found: {clim_file}"
    clm_time = np.linspace(1, 365.25, 12)

    time_ref_list = list(clm_time)
    info(f"Loaded climatology file: {os.path.basename(clim_file)}")
    info(f"time dimension = {len(time_ref_list)} (should be 12 months)")

    # OGCM metadata (reference = SSH/zeta)
    ogcm = tl.load_ogcm_metadata(clim_file, cfg.bio_ogcm_var_name)

# --- [04] Prepare weights ---
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.weight_file}")
        plus(f"Weight file created: {cfg.weight_file}")
    else:
        info(f"Use existing wght file {cfg.weight_file}")

    with Dataset(cfg.weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S   = nc.variables["S"][:]

# Prepare bry data containers
with step("[04b] Allocate bry buffers"):
    bry_data = tl.make_all_bry_data_shapes(
        ["NO3","PO4","NH4","TIC","oxygen","chlo","phyt","zoop","alkalinity","SdeC","LdeC","RdeC","SdeN","RdeN"],
        len(time_ref_list),
        grd,
        cfg.vertical.layer_n,
    )
    bry_time = []

# --- [05] Main processing loop (12 climatology months) ---
nc_clim = tl.MaskedNetCDF(Dataset(cfg.bio_ogcm_name, mode="r", maskandscale=True))

for i, t in enumerate(time_ref_list):
    tag = f"month={i+1:02d}"

    with step("[07] Load OGCM fields", ts=tag):
        NO3  = nc_clim.get(cfg.bio_ogcm_var_name['NO3'], i, slice(None), idy, idx)
        phyt = nc_clim.get(cfg.bio_ogcm_var_name['phyt'], i, slice(None), idy, idx)
        chlo = nc_clim.get(cfg.bio_ogcm_var_name['chl'], i, slice(None), idy, idx)
        PO4 = nc_clim.get(cfg.bio_ogcm_var_name['PO4'], i, slice(None), idy, idx)
        oxygen = nc_clim.get(cfg.bio_ogcm_var_name['oxygen'], i, slice(None), idy, idx)
        zoop = phyt*0.3

        alkalinity = np.ones_like(zoop) * 2350
        TIC = np.ones_like(zoop) * 2100
        NH4 = np.ones_like(zoop) * 0.01

        SdeC = np.ones_like(zoop) * 0.04
        LdeC = np.ones_like(zoop) * 0.04
        RdeC = np.ones_like(zoop) * 0.04
        SdeN = np.ones_like(zoop) * 0.04
        RdeN = np.ones_like(zoop) * 0.04

        field = tl.ConfigObject(NO3=NO3, phyt=phyt, chlo=chlo, PO4=PO4, oxygen=oxygen, zoop=zoop,alkalinity=alkalinity, TIC=TIC, NH4=NH4, SdeC=SdeC, LdeC=LdeC,RdeC=RdeC, SdeN=SdeN, RdeN=RdeN)
        
    # --- 나머지 [08] Remap → [11] z→σ 저장 로직은 기존 코드 그대로 ---
    #    field에 대해 remap/flood/mask/rotate/sigma transform 동일하게 실행

        # [08] Remap (weights)
        with step("[08] Remap (weights)"):
            for var in vars(field):
                var_src = getattr(field, var)
                with capture_warnings(tag):
                    remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
                setattr(field, var, remapped)

        # [09] Flood H/V
        with step("[09] Flood H/V", ts=f"None"):
            for var in vars(field):
                val = getattr(field, var)
                val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_bry)
                setattr(field, var, val_flooded)
            for var in vars(field):
                val = getattr(field, var)
                if val.ndim == 2:
                    continue
                with capture_warnings(tag):
                    val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask), spval=-1e10)
                setattr(field, var, val_flooded)

        # [10] Mask 
        with step("[10] Mask & rotate", ts=f"None"):
            for var in vars(field):
                arr = getattr(field, var)
                if arr.ndim == 2:
                    arr[grd.mask == 0] = 0.0
                else:
                    arr[..., grd.mask == 0] = 0.0
                setattr(field, var, arr)

        # [11] z→σ & save bry
        with step("[11] z→σ & save bry", ts=f"None"):
            directions = ['north', 'south', 'east', 'west']

            # HYCOM depth padding
            Z = np.zeros(len(ogcm.depth) + 2)
            Z[0] = 100
            Z[1:-1] = -np.abs(ogcm.depth)
            Z[-1] = -100000
            Z_flipped = np.flipud(Z)

            # sigma related parameters
            zargs = (
                cfg.vertical.vtransform,
                cfg.vertical.vstretching,
                cfg.vertical.theta_s,
                cfg.vertical.theta_b,
                cfg.vertical.tcline,
                cfg.vertical.layer_n,
            )

            zr_3d = tl.zlevs(*zargs, 1, grd.topo, np.zeros_like(grd.topo))    # (Ns, Mp, Lp)
            zw_3d = tl.zlevs(*zargs, 5, grd.topo, np.zeros_like(grd.topo))    # (Ns+1, Mp, Lp)
            dz_3d = zw_3d[1:, :, :] - zw_3d[:-1, :, :]


            for direction in directions:
                zr = tl.extract_bry(zr_3d, direction)
                dzr = tl.extract_bry(dz_3d, direction)

                for var, zgrid in zip(["NO3","PO4","NH4","TIC","oxygen","chlo","phyt","zoop","alkalinity","SdeC","LdeC","RdeC","SdeN","RdeN"], [zr, zr, zr, zr,zr,zr,zr,zr,zr,zr,zr,zr,zr,zr]):
                    val = tl.extract_bry(getattr(field, var), direction)        # (Nz, Lp)
                    val_pad = np.vstack((val[0:1], val, val[-1:]))              # (Nz+2, Lp)
                    val_flip = np.flip(val_pad, axis=0)                         # (Nz+2)
                   
                    with capture_warnings(tag):
                        sigma_interped = tl.ztosigma_1d_numba(val_flip, zgrid, Z_flipped)  # (Ns, Lp)
                    barotropic = None

                    if barotropic is not None:
                        with capture_warnings(tag):
                            sigma_interped, _, barotropic_corrected, _ = tl.conserve_and_recompute_barotropic(
                                sigma_interped, sigma_interped, barotropic, barotropic, dz_bar, dz_bar
                            )

                    setattr(field, f"{var}_{direction}", sigma_interped)

                # save to buffers
                for varname in bry_data:
                    val_d = getattr(field, f"{varname}_{direction}")
                    bry_data[varname][direction][i, ...] = val_d
                   
                # time save

        # per-timestep done line

# [12] Write all remapped variables to bry.nc
with step("[12] Write variables", out=ellipsis(cfg.bryname)):
    with Dataset(cfg.bryname, 'a') as nc:
        for varname in bry_data:
            for direction in bry_data[varname]:
                var_fullname = f"{varname}_{direction}"
                nc[var_fullname][:] = bry_data[varname][direction]

bar("Summary")
print(f"Total elapsed: {time.time() - start1:.3f}s")

