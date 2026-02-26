
# add_ini_fennel_clm.py (refactored)
#
# - bio variable list(ini 기준)을 bio_vars.yml에서 읽어서 자동으로 처리
# - cfg.bio_ogcm_var_name 매핑에 있는 것만 OGCM에서 읽고, 없는 건 0 유지(또는 아래에서 파생/상수로 채움)
# - remap / flood / z→σ / write 단계는 “field에 들어있는 변수”만 대상으로 수행
#
# 원본: add_ini_fennel_clm.py :contentReference[oaicite:0]{index=0}

# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
from netCDF4 import Dataset, date2num
import time

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..", "libs")))

import create_I as cn
import utils as tl

# YAML bio schema loader: libs/create.py 쪽 함수 재사용 (없으면 로컬 fallback)
try:
    import create as cr  # libs/create.py
except Exception:  # 아주 예외적인 경우만
    cr = None
    import yaml

    def _load_bio_yaml_fallback(path: str) -> dict:
        try:
            with open(path, "r") as f:
                return yaml.safe_load(f) or {}
        except Exception:
            return {}

    def _get_bio_defs_fallback(bio_db: dict, bio_model, target: str) -> dict:
        if bio_model is None:
            return {}
        key = str(bio_model).lower()
        models = (bio_db or {}).get("bio_models", {})
        return (models.get(key, {}).get(target, {}) or {})

# logging
from log_utils2 import configure, step, capture_warnings, info, note, bar, warn_line, ellipsis

configure(width=80, show_sections=False, color_mode="auto")


def _get_bio_names_from_yaml(cfg) -> list[str]:
    """
    bio_vars.yml에서 ini용 bio 변수명 리스트를 가져온다.
    - 파일 없거나 YAML 깨지면: 빈 리스트(=bio 스킵)
    """
    bio_yaml = getattr(cfg, "bio_yaml", "bio_vars.yml")
    bio_model = getattr(cfg, "bio_model_type", None)

    if cr is not None:
        bio_db = cr.load_bio_yaml(bio_yaml)
        defs = cr.get_bio_defs(bio_db, bio_model, target="ini")
    else:
        bio_db = _load_bio_yaml_fallback(bio_yaml)
        defs = _get_bio_defs_fallback(bio_db, bio_model, target="ini")

    return list((defs or {}).keys())


def _first_3d_template(field_obj):
    """
    상수필드(알칼리니티/TIC/NH4 등) 만들 때 shape 템플릿 잡으려고,
    field 안에서 3D (z, y, x) 배열 하나 찾는다.
    """
    for name in vars(field_obj):
        arr = getattr(field_obj, name)
        if isinstance(arr, np.ndarray) and arr.ndim == 3:
            return arr
    return None


# --- main --------------------------------------------------------------------
start0 = time.time()
bar("Initial Condition Build (ini)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg = tl.parse_config("./config.yaml")
    grd = tl.load_roms_grid(cfg.grdname)
    ogcm = tl.load_ogcm_metadata(cfg.bio_ogcm_name, cfg.bio_ogcm_var_name)
    info(f"grid={cfg.grdname}")
    info(f"ini_out={cfg.ininame}")
    info(f"wght_out={cfg.bio_weight_file}")
    info(f"ogcm={cfg.bio_ogcm_name}")

with step("[02] Load climatology file"):
    clim_file = cfg.bio_ogcm_name
    assert os.path.exists(clim_file), f"Climatology file not found: {clim_file}"
    relative_time = date2num(cfg.initdate, cfg.time_ref)
    info(f"Loaded climatology file: {os.path.basename(clim_file)}")

# [03] Create initial NetCDF file (원하면 다시 켜)
# with step("[03] Create initial NetCDF"):
#     status = cn.create_ini__(cfg, grd, relative_time, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
#     if status:
#         raise RuntimeError(f"Failed creating file {cfg.ininame}")

# [04] Crop OGCM domain and prepare remap weights
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.bio_weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.bio_weight_file}")
    else:
        info(f"Use existing wght file {cfg.bio_weight_file}")

# --- bio variable list from YAML (ini target) ---
with step("[05] Resolve bio variable list (YAML)"):
    bio_names = _get_bio_names_from_yaml(cfg)
    if not bio_names:
        note("bio_vars.yml에서 ini bio 변수를 못 읽어서(또는 bio_model 없음) 물리만 진행할 수 있음")
    else:
        info(f"bio_model={cfg.bio_model_type} | ini bio vars={len(bio_names)}")

nc_clim = tl.MaskedNetCDF(Dataset(cfg.bio_ogcm_name, mode="r", maskandscale=True))

# [06] Load OGCM fields (YAML 기반)
with step("[06] Load OGCM fields", ts="Load OGCM"):
    field_dict = {}

    # 6-1) OGCM에서 직접 읽을 수 있는 것만 읽는다 (cfg.bio_ogcm_var_name에 매핑이 있는 경우)
    # 매핑 키는 "ini 변수명"과 동일해야 한다.
    for vname in bio_names:
        src_name = cfg.bio_ogcm_var_name.get(vname, None)
        if src_name is None:
            warn_line(f"OGCM mapping missing: cfg.bio_ogcm_var_name['{vname}'] (keep 0)")
            continue
        field_dict[vname] = nc_clim.get(src_name, 0, slice(None), idy, idx)

    # 6-2) 파생/상수로 채우고 싶은 변수들(있을 때만)
    # - 여기 규칙은 “있으면 채움, 없으면 스킵”
    # - OGCM에서 못 읽은 변수는 최종적으로 0 유지될 수 있음(ini 껍데기 정책상 OK)

    # zooplankton = phytoplankton * 0.3 (둘 다 YAML에 있을 때만)
    if ("zooplankton" in bio_names) and ("phytoplankton" in field_dict):
        field_dict["zooplankton"] = field_dict["phytoplankton"] * 0.3

    # 상수 필드는 템플릿 shape이 필요함
    field = tl.ConfigObject(**field_dict)
    tmpl3 = _first_3d_template(field)

    def _fill_const(name: str, value: float):
        if name not in bio_names:
            return
        if name in field_dict:
            return
        if tmpl3 is None:
            warn_line(f"cannot create constant field '{name}': no 3D template loaded")
            return
        field_dict[name] = np.ones_like(tmpl3) * value

    _fill_const("alkalinity", 2350.0)
    _fill_const("TIC", 2100.0)
    _fill_const("NH4", 0.01)

    _fill_const("SdetritusC", 0.04)
    _fill_const("LdetritusC", 0.04)
    _fill_const("RdetritusC", 0.04)
    _fill_const("SdetritusN", 0.04)
    _fill_const("RdetritusN", 0.04)

    # field 갱신
    field = tl.ConfigObject(**field_dict)

# [07] Load and apply remap weights to all fields
with step("[07] Remap (weights)"):
    with Dataset(cfg.bio_weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S = nc.variables["S"][:]

    for varname in vars(field):
        var_src = getattr(field, varname)
        with capture_warnings(tag="remap"):
            remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
        setattr(field, varname, remapped)

note(f"Flood method: {cfg.flood_method_for_ini}")

# [08] Horizontal flood (all fields)
with step("[08] Flood: horizontal"):
    for var in vars(field):
        val = getattr(field, var)
        val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_ini)
        setattr(field, var, val_flooded)

# [09] Vertical flood
with step("[09] Flood: vertical"):
    for var in vars(field):
        val = getattr(field, var)
        with capture_warnings(tag="vflood"):
            val_flooded = tl.flood_vertical_vectorized(val, grd.mask_rho, spval=-1e10)
        setattr(field, var, val_flooded)

# [10] Mask land to 0 (원본 라벨 유지)
with step("[10] Mask & rotate", ts="None"):
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

    Z = np.zeros(len(ogcm.depth) + 2)
    Z[0] = 100
    Z[1:-1] = -np.abs(ogcm.depth)
    Z[-1] = -100000
    Zf = np.flipud(Z)

    # field에 있는 3D 변수만 z→σ 수행
    for var in vars(field):
        val = getattr(field, var)
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
    with Dataset(cfg.ininame, "a") as nc:
        for varname in vars(field):
            if varname not in nc.variables:
                warn_line(f"skip write (not in ini): {varname}")
                continue
            nc[varname][0] = getattr(field, varname)

bar("Summary")
