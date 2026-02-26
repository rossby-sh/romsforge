
# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
from netCDF4 import Dataset, date2num
import time

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl
from io_utils import load_bio_yaml, get_bio_vars, get_bio_rules, apply_bio_rules

# logging
from log_utils2 import configure, step, capture_warnings, info, note, plus, bar, warn_line, done, ellipsis
configure(width=80, show_sections=False, color_mode='auto')

# --- [01] Helpers ----------------------------------------------------------------
def pick_daily_ogcm_file(cfg):
    """
    cfg.bio_ogcm_name: daily directory (cmems_bio_YYYYMMDD.nc files)
    return: full path to cmems_bio_YYYYMMDD.nc for cfg.initdate
    """
    daily_dir = cfg.bio_ogcm_name
    if not os.path.isdir(daily_dir):
        raise NotADirectoryError(f"cfg.bio_ogcm_name must be a directory: {daily_dir}")

    ymd = f"{cfg.initdate:%Y%m%d}"
    ogcm_file = os.path.join(daily_dir, f"cmems_bio_{ymd}.nc")
    if not os.path.exists(ogcm_file):
        raise FileNotFoundError(ogcm_file)
    return ogcm_file


# --- main --------------------------------------------------------------------
start0 = time.time()
bar("Initial Condition Build (ini)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg  = tl.parse_config("./config.yaml")
    grd  = tl.load_roms_grid(cfg.grdname)

    # cfg.bio_ogcm_name is a DAILY DIRECTORY
    ogcm_file = pick_daily_ogcm_file(cfg)

    # load metadata from the picked daily file
    ogcm = tl.load_ogcm_metadata(ogcm_file, cfg.bio_ogcm_var_name)

    # IMPORTANT: cfg.bio_ogcm_var_name is ConfigObject -> use dict view for membership/get
    bio_map = cfg.bio_ogcm_var_name.__dict__

    info(f"grid={cfg.grdname}")
    info(f"ini_out={cfg.ininame}")
    info(f"wght_out={cfg.bio_weight_file}")
    info(f"ogcm(daily)={ogcm_file}")
    info(f"bio_model={cfg.bio_model_type}")

with step("[02] Load daily ogcm file"):
    assert os.path.exists(ogcm_file), f"Daily OGCM file not found: {ogcm_file}"
    relative_time = date2num(cfg.initdate, cfg.time_ref)
    info(f"Loaded daily ogcm file: {os.path.basename(ogcm_file)}")

# [03] Create initial NetCDF file
#with step("[03] Create initial NetCDF"):
#    status = cn.create_ini__(cfg, grd, relative_time, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
#    if status:
#        raise RuntimeError(f"Failed creating file {cfg.ininame}")

# [03b] Load bio_vars.yml (ini vars + rules)
with step("[03b] Load bio vars (YAML)"):
    bio_yaml = getattr(cfg, "bio_yaml", "bio_vars.yml")
    bio_db = load_bio_yaml(bio_yaml)

    bio_model = str(cfg.bio_model_type).lower()
    ini_vars = get_bio_vars(bio_db, bio_model, target="ini")
    bio_rules = get_bio_rules(bio_db, bio_model,target="ini")

    info(f"bio_yaml={bio_yaml}")
    info(f"ini_vars={ini_vars}")

# [04] Crop OGCM domain and prepare remap weights
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.bio_weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.bio_weight_file}")
        plus(f"Weight file created: {cfg.bio_weight_file}")
    else:
        info(f"Use existing wght file {cfg.bio_weight_file}")

# [05] Load OGCM fields (yaml-driven) + rules
nc_clim = tl.MaskedNetCDF(Dataset(ogcm_file, mode="r", maskandscale=True))
with step("[05] Load OGCM fields", ts="Load OGCM"):
    field_dict = {}

    # (1) OGCM에서 직접 읽을 수 있는 것만 로딩
    for var in ini_vars:
        src = bio_map.get(var, None)
        if src is None:
            continue
        field_dict[var] = nc_clim.get(src, 0, slice(None), idy, idx)

    # >>> FIXME >>>
    print("!!! SHJO: FIX ME !!!")
    field_dict["phytoplankton"]=field_dict["phytoplankton"]/(0.02*6.625*12) 
    # <<< FIXME <<<

    # (2) rules로 파생/상수 채우기
    apply_bio_rules(field_dict, bio_rules)

    # (3) 누락 체크
    missing = [v for v in ini_vars if v not in field_dict]
    if missing:
        raise KeyError(f"ini missing vars (no mapping & no rule): {missing}")

    field = tl.ConfigObject(**field_dict)

# [06] Load and apply remap weights to all fields
with step("[06] Remap (weights)"):
    with Dataset(cfg.bio_weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S   = nc.variables["S"][:]

    for varname in ini_vars:
        var_src = getattr(field, varname)
        with capture_warnings(tag="remap"):
            remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
        setattr(field, varname, remapped)



note(f"Flood method: {cfg.flood_method_for_ini}")

# [07] Horizontal flood (all fields)
with step("[07] Flood: horizontal"):
    for var in ini_vars:
        val = getattr(field, var)
        val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_ini)
        setattr(field, var, val_flooded)

# [08] Vertical flood
with step("[08] Flood: vertical"):
    for var in ini_vars:
        val = getattr(field, var)
        with capture_warnings(tag="vflood"):
            val_flooded = tl.flood_vertical_vectorized(val, grd.mask_rho, spval=-1e10)
        setattr(field, var, val_flooded)

# [09] Mask land to 0
with step("[10] Mask & rotate", ts=f"None"):
    for var in ini_vars:
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

    Z = np.zeros(len(ogcm.depth) + 2)
    Z[0] = 100
    Z[1:-1] = -np.abs(ogcm.depth)
    Z[-1] = -100000
    Zf = np.flipud(Z)

    # 생물 변수는 전부 rho-point 가정 → zr 고정
    for var in ini_vars:
        val = getattr(field, var)

        # 3D만 z→σ 수행 (Nz, Mp, Lp)
        if not isinstance(val, np.ndarray) or val.ndim != 3:
            continue

        padded = np.vstack((val[0:1], val, val[-1:]))
        flipped = np.flip(padded, axis=0)
        with capture_warnings(tag="z2sigma"):
            remapped = tl.ztosigma_numba(flipped, zr, Zf)
        setattr(field, var, remapped)

bar("Summary")
print(f"Total elapsed: {time.time() - start0:.3f}s")

# [12] Write all remapped variables to ini.nc
with step("[12] Write variables", out=ellipsis(cfg.ininame)):
    with Dataset(cfg.ininame, 'a') as nc:
        # 시간 맞추고 싶으면 켜
        # if "ocean_time" in nc.variables:
        #     nc["ocean_time"][0] = relative_time

        for varname in ini_vars:
            if varname not in nc.variables:
                warn_line(f"skip write (not in ini): {varname}")
                continue
            nc[varname][0] = getattr(field, varname)

bar("Summary")
