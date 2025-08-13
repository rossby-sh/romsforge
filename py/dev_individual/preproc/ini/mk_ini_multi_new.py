# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
from netCDF4 import Dataset, num2date, date2num
import time
from collections import defaultdict

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_I as cn
import utils as tl
from io_utils import collect_time_info  # filelist 기반

# logging
from log_utils2 import configure, step, capture_warnings, info, note, plus, warn_line, done, bar
configure(width=80, show_sections=False, color_mode='auto')

# -------- helpers -------------------------------------------------------------
def resolve_filelist(path_or_pattern: str):
    """Accepts directory, exact file, or glob pattern. Returns a sorted list."""
    p = path_or_pattern
    # explicit file
    if os.path.isfile(p):
        return [p]
    # directory -> *.nc
    if os.path.isdir(p):
        return sorted(glob.glob(os.path.join(p, "*.nc")))
    # glob pattern
    lst = sorted(glob.glob(p))
    return lst

# --- main --------------------------------------------------------------------
start0 = time.time()
bar("Initial Condition Build (ini, multi-source)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg  = tl.parse_config("./config_multi.yaml")
    grd  = tl.load_roms_grid(cfg.grdname)
    ref_key = getattr(cfg, "reference_var", "zeta")
    info(f"grid={cfg.grdname}")
    info(f"ini_out={cfg.ininame}")
    info(f"wght_out={cfg.weight_file}")
    note(f"reference_var={ref_key}")

# 멀티 소스: 변수별 파일/시간 수집
filename_dict = {}     # var -> matched file list (ordered by time)
tinfo_dict = {}        # var -> tinfo list
time_ref_list = None   # 기준 시간축

# [02] Index OGCM files & match init time per variable
with step("[02] Index OGCM files & match init time"):
    for varname, meta in cfg.ogcm_inputs.to_dict().items():
        path = meta["path"]
        var  = meta["varname"]

        filelist = resolve_filelist(path)
        assert len(filelist) > 0, f"No files found for variable: {varname} (path='{path}')"

        # initdate 하나만 매칭 (start=end=initdate)
        tinfo = collect_time_info(
            filelist,
            cfg.ogcm_var_name["time"],
            (str(cfg.initdate), str(cfg.initdate)),
        )
        assert len(tinfo) >= 1, f"No time match at initdate for variable: {varname}"

        tinfo_dict[varname] = tinfo
        filename_dict[varname] = [ti.filename for ti in tinfo]
        info(f"{varname}: files_matched={len(tinfo)} first={os.path.basename(filename_dict[varname][0])}")

    # 모든 변수의 시간 일치성 검사
    for varname, tlist in tinfo_dict.items():
        dt_list = [ti.datetime for ti in tlist]
        if time_ref_list is None:
            time_ref_list = dt_list
        else:
            if time_ref_list != dt_list:
                raise ValueError(f"Time mismatch in variable {varname}")
    # 단일 initdate가 이상적
    assert len(time_ref_list) == 1, f"Expected a single init time, got {len(time_ref_list)}"
    note(f"All variables aligned at init: {time_ref_list[0].isoformat()}")

# [02b] Load OGCM metadata from reference stream
with step("[02b] Load OGCM metadata (ref)"):
    if ref_key not in filename_dict:
        raise KeyError(f"reference_var='{ref_key}' not found in ogcm_inputs")
    ref_file = filename_dict[ref_key][0]
    ogcm = tl.load_ogcm_metadata(ref_file, cfg.ogcm_var_name)
    ref_tinfo = tinfo_dict[ref_key][0]
    relative_time = tl.compute_relative_time(ref_tinfo.raw_value, ogcm.time_unit, cfg.time_ref)
    info(f"ref_file={ref_file}")

# [03] Create initial NetCDF file
with step("[03] Create initial NetCDF"):
    status = cn.create_ini__(cfg, grd, relative_time, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
    if status:
        raise RuntimeError(f"Failed creating file {cfg.ininame}")
    plus(f"Created file: {cfg.ininame}")

# [04] Crop OGCM domain and prepare remap weights
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.weight_file}")
        plus(f"Weight file created: {cfg.weight_file}")
    else:
        info(f"Use existing wght file {cfg.weight_file}")

# [05] Load OGCM fields at init (multi inputs)
with step("[05] Load OGCM fields (@init)"):
    sel = {v: tinfo_dict[v][0] for v in tinfo_dict.keys()}
    tag = f"ts={time_ref_list[0].strftime('%Y-%m-%d %H')}"
    with Dataset(sel['zeta'].filename, maskandscale=True) as ncz, \
         Dataset(sel['temp'].filename, maskandscale=True) as nct, \
         Dataset(sel['salt'].filename, maskandscale=True) as ncs, \
         Dataset(sel['u'].filename,    maskandscale=True) as ncu, \
         Dataset(sel['v'].filename,    maskandscale=True) as ncv:
        nz = tl.MaskedNetCDF(ncz); nt = tl.MaskedNetCDF(nct)
        ns = tl.MaskedNetCDF(ncs); nu = tl.MaskedNetCDF(ncu); nv = tl.MaskedNetCDF(ncv)
        iz = sel['zeta'].index; it = sel['temp'].index; isal = sel['salt'].index; iu = sel['u'].index; iv = sel['v'].index

        zeta = nz.get(cfg.ogcm_inputs['zeta']['varname'], iz,               idy, idx)
        temp = nt.get(cfg.ogcm_inputs['temp']['varname'], it, slice(None),  idy, idx)
        salt = ns.get(cfg.ogcm_inputs['salt']['varname'], isal, slice(None), idy, idx)
        u    = nu.get(cfg.ogcm_inputs['u']['varname'],    iu, slice(None),  idy, idx)
        v    = nv.get(cfg.ogcm_inputs['v']['varname'],    iv, slice(None),  idy, idx)
        with capture_warnings(tag=tag):
            ubar = tl.depth_average(u, ogcm.depth)
            vbar = tl.depth_average(v, ogcm.depth)

    field = tl.ConfigObject(zeta=zeta, ubar=ubar, vbar=vbar, temp=temp, salt=salt, u=u, v=v)

# [06] Remap (weights)
with step("[06] Remap (weights)"):
    with Dataset(cfg.weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S   = nc.variables["S"][:]
    for varname in vars(field):
        var_src = getattr(field, varname)
        with capture_warnings(tag="remap"):
            remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
        setattr(field, varname, remapped)

note(f"Flood method: {cfg.flood_method_for_ini}")

# [07] Flood: horizontal
with step("[07] Flood: horizontal"):
    for var in ['temp', 'salt', 'u', 'v','zeta','ubar','vbar']:
        val = getattr(field, var)
        val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_ini)
        setattr(field, var, val_flooded)

# [08] Flood: vertical
with step("[08] Flood: vertical"):
    for var in ['temp', 'salt', 'u', 'v']:
        val = getattr(field, var)
        with capture_warnings(tag="vflood"):
            val_flooded = tl.flood_vertical_vectorized(val, grd.mask, spval=-1e10)
        setattr(field, var, val_flooded)

# [09] Mask & clean
with step("[09] Mask & clean"):
    for varname in ['zeta', 'ubar', 'vbar']:
        var = getattr(field, varname)
        var[grd.mask == 0] = 0.0
        setattr(field, varname, var)
    for varname in ['temp', 'salt', 'u', 'v']:
        var = getattr(field, varname)
        var[:, grd.mask == 0] = 0.0
        setattr(field, varname, var)

# [10] Rotate & stagger (u,v)
with step("[10] Rotate & stagger (u,v)"):
    u_rot, v_rot       = tl.rotate_vector_euler(field.u,    field.v,    grd.angle, to_geo=False)
    ubar_rot, vbar_rot = tl.rotate_vector_euler(field.ubar, field.vbar, grd.angle, to_geo=False)
    setattr(field, 'u',    tl.rho2uv(u_rot,'u'))
    setattr(field, 'v',    tl.rho2uv(v_rot,'v'))
    setattr(field, 'ubar', tl.rho2uv(ubar_rot,'u'))
    setattr(field, 'vbar', tl.rho2uv(vbar_rot,'v'))

# [11] z→σ interpolation
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
    Z[0] = 100; Z[1:-1] = -ogcm.depth; Z[-1] = -100000
    Zf = np.flipud(Z)

    for var, zgrid in zip(['temp', 'salt', 'u', 'v'], [zr, zr, zu, zv]):
        val = getattr(field, var)
        padded  = np.vstack((val[0:1], val, val[-1:]))
        flipped = np.flip(padded, axis=0)
        with capture_warnings(tag="z2sigma"):
            remapped = tl.ztosigma_numba(flipped, zgrid, Zf)
        setattr(field, var, remapped)

# [12] Conserve volume & fix barotropic
with step("[12] Conserve volume & fix barotropic"):
    zw  = tl.zlevs(*zlevs_args, 5, grd.topo, field.zeta)
    dzr = zw[1:, :, :] - zw[:-1, :, :]
    dzu, dzv = tl.rho2uv(dzr,'u'), tl.rho2uv(dzr,'v')

    with capture_warnings(tag="baro"):
        u_consv, v_consv, ubar_new, vbar_new = tl.conserve_and_recompute_barotropic(
            field.u, field.v, field.ubar, field.vbar, dzu, dzv
        )
    setattr(field, 'u', u_consv)
    setattr(field, 'v', v_consv)
    setattr(field, 'ubar', ubar_new)
    setattr(field, 'vbar', vbar_new)

# [13] Write variables
with step("[13] Write variables", out=cfg.ininame):
    with Dataset(cfg.ininame, mode='a') as nc:
        for var in ['zeta', 'temp', 'salt', 'u', 'v', 'ubar', 'vbar']:
            nc[var][0] = getattr(field, var)

bar("Summary")
print(f"Total elapsed: {time.time() - start0:.3f}s")

