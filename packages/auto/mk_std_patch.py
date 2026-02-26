# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
from netCDF4 import Dataset, num2date, date2num
import time

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create as cn
import utils as tl
from io_utils import collect_time_info_legacy

# logging
from log_utils2 import configure, step, capture_warnings, info, note, plus, warn_line, done, bar
configure(width=80, show_sections=False, color_mode='auto')

# --- main --------------------------------------------------------------------
start0 = time.time()
bar("Initial Condition Build (ini)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg  = tl.parse_config("./config_da.yaml")
    grd  = tl.load_roms_grid(cfg.grdname)
    ogcm = tl.load_ogcm_metadata(cfg.ogcm_name, cfg.ogcm_var_name)
    info(f"grid={cfg.grdname}")
    info(f"wght_out={cfg.weight_file}")
    info(f"ogcm={cfg.ogcm_name}")

# [02] Time index matching & relative time
#with step("[02] Time index matching & relative time"):
#    tinfo = collect_time_info_legacy(cfg.ogcm_name, cfg.ogcm_var_name['time'], str(cfg.initdate))
#    _, idt, _, _ = tinfo[0]
#    relative_time = tl.compute_relative_time(ogcm.time[idt], ogcm.time_unit, cfg.time_ref)

# [03] Create STD NetCDF files (created per-month in loop)


# [04] Crop OGCM domain and prepare remap weights
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.weight_file}")
    else:
        info(f"Use existing wght file {cfg.weight_file}")


tvec = np.arange(15, 366, 365.25/12)   # 12개 나옴
for idt in range(12):

    # 월별 시간: 365.25를 12등분 (days)
    relative_time = float(tvec[idt])*86400

    month_str_num = f"{idt+1:02d}"
    month_str_name = ["jan","feb","mar","apr","may","jun",
                      "jul","aug","sep","oct","nov","dec"][idt]

    if cfg.month_num:
        month_tag = month_str_num
    else:
        month_tag = month_str_name

    outname = os.path.join(
        cfg.std_dir,
        f"{cfg.std_file}_{month_tag}.nc"
    )


    # [03] Create STD NetCDF (single month file)
    with step("[03] Create std NetCDF", out=outname):
        status = cn.create_std(cfg, grd, relative_time, outname, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
        if status:
            raise RuntimeError(f"Failed creating file {outname}")

# [05] Load OGCM raw fields (single time)

    with step("[05] Load OGCM fields"):
        t0 = time.time()
        with Dataset(cfg.ogcm_name, maskandscale=True) as nc_raw:
            nc = tl.MaskedNetCDF(nc_raw)
            zeta = nc.get(cfg.ogcm_var_name['zeta'],        idt,              idy, idx)
            temp = nc.get(cfg.ogcm_var_name['temperature'], idt, slice(None), idy, idx)
            salt = nc.get(cfg.ogcm_var_name['salinity'],    idt, slice(None), idy, idx)
            u    = nc.get(cfg.ogcm_var_name['u'],           idt, slice(None), idy, idx)
            v    = nc.get(cfg.ogcm_var_name['v'],           idt, slice(None), idy, idx)
            with capture_warnings(tag=f"ts={getattr(cfg, 'initdate', 'NA')}"):
                ubar = tl.depth_average(u, ogcm.depth)
                vbar = tl.depth_average(v, ogcm.depth)
        field = tl.ConfigObject(zeta=zeta, ubar=ubar, vbar=vbar, temp=temp, salt=salt, u=u, v=v)
#    done("load", time.time()-t0)

# [06] Load and apply remap weights to all fields
    with step("[06] Remap (weights)"):
        with Dataset(cfg.weight_file) as nc:
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



    for var in ['temp','u','v','salt']:
        val = getattr(field, var)
        if val.ndim == 3:
            n0 = np.count_nonzero(~np.isnan(val[0]))
            nmid = np.count_nonzero(~np.isnan(val[val.shape[0]//2]))
            nmin = min(np.count_nonzero(~np.isnan(val[k])) for k in range(val.shape[0]))
            print(var, "valid:", "surf", n0, "mid", nmid, "min_over_k", nmin)
        else:
            print(var, "valid:", np.count_nonzero(~np.isnan(val)))

    note(f"Flood method: {cfg.flood_method_for_ini}")
# [07] Horizontal flood (all fields)
    with step("[07] Flood: horizontal"):
        t0 = time.time()
        for var in ['temp', 'salt', 'u', 'v','zeta','ubar','vbar']:
            val = getattr(field, var)
            val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_ini)
            setattr(field, var, val_flooded)
#    done("flood_h", time.time()-t0)

# [08] Vertical flood
    with step("[08] Flood: vertical"):
        t0 = time.time()
        for var in ['temp', 'salt', 'u', 'v']:
            val = getattr(field, var)
            with capture_warnings(tag="vflood"):
                val_flooded = tl.flood_vertical_vectorized(val, grd.mask_rho, spval=-1e10)
                # Alternative:
                # val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask), spval=-1e10)
            setattr(field, var, val_flooded)
#    done("flood_v", time.time()-t0)

# [09] Mask land to 0
    with step("[09] Mask & clean"):
        for varname in ['zeta', 'ubar', 'vbar']:
            var = getattr(field, varname)
            var[grd.mask_rho == 0] = 0.0
            setattr(field, varname, var)
        for varname in ['temp', 'salt', 'u', 'v']:
            var = getattr(field, varname)
            var[:, grd.mask_rho == 0] = 0.0
            setattr(field, varname, var)

# [10] Rotate vectors to model grid and convert to u/v points
    with step("[10] Rotate & stagger (u,v)"):
        u_rot, v_rot       = tl.rotate_vector_euler(field.u,    field.v,    grd.angle, to_geo=False)
        ubar_rot, vbar_rot = tl.rotate_vector_euler(field.ubar, field.vbar, grd.angle, to_geo=False)
        setattr(field, 'u',    tl.rho2uv(u_rot,'u'))
        setattr(field, 'v',    tl.rho2uv(v_rot,'v'))
        setattr(field, 'ubar', tl.rho2uv(ubar_rot,'u'))
        setattr(field, 'vbar', tl.rho2uv(vbar_rot,'v'))

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
        zr = tl.zlevs(*zlevs_args, 1, grd.topo, field.zeta)
        zu, zv = tl.rho2uv(zr,'u'), tl.rho2uv(zr,'v')

        Z = np.zeros(len(ogcm.depth)+2)
        Z[0] = 100; Z[1:-1] = -np.abs(ogcm.depth); Z[-1] = -100000
        Zf = np.flipud(Z)

        for var, zgrid in zip(['temp', 'salt', 'u', 'v'], [zr, zr, zu, zv]):
            val = getattr(field, var)
            padded = np.vstack((val[0:1], val, val[-1:]))
            flipped = np.flip(padded, axis=0)
            with capture_warnings(tag="z2sigma"):
                remapped = tl.ztosigma_numba(flipped, zgrid, Zf)
            setattr(field, var, remapped)

# [13] Write all variables to ini.nc
    with step("[13] Write variables", out=outname):
        with Dataset(outname, mode='a') as nc:
            for var in ['zeta', 'temp', 'salt', 'u', 'v', 'ubar', 'vbar']:
                nc[var][0] = getattr(field, var) * cfg.scale_factor

    bar("Summary")
    print(f"Total elapsed: {time.time() - start0:.3f}s")

