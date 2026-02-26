
# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
from netCDF4 import Dataset
import time

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl
from io_utils import load_bio_yaml, get_bio_vars, get_bio_rules, apply_bio_rules

# logging
from log_utils2 import configure, step, capture_warnings, info, note, plus, bar, warn_line, done, ellipsis
configure(width=80, show_sections=False, color_mode='auto')

# --- [01] Helpers ----------------------------------------------------------------
def month_tag(idt: int, cfg):
    """Return month tag for output filename."""
    months = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]
    if getattr(cfg, "month_num", True):
        return f"{idt+1:02d}"
    return months[idt]


# --- main --------------------------------------------------------------------
start0 = time.time()
bar("STD Build (bio)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg  = tl.parse_config("./config_da.yaml")
    grd  = tl.load_roms_grid(cfg.grdname)

    # monthly/std OGCM file (single file with time=12)
    ogcm_file = cfg.bio_ogcm_name
    if not os.path.exists(ogcm_file):
        raise FileNotFoundError(ogcm_file)

    ogcm = tl.load_ogcm_metadata(ogcm_file, cfg.bio_ogcm_var_name)

    # IMPORTANT: cfg.bio_ogcm_var_name is ConfigObject -> use dict view for membership/get
    bio_map = cfg.bio_ogcm_var_name.__dict__

    info(f"grid={cfg.grdname}")
    info(f"std_dir={cfg.std_dir}")
    info(f"std_file={cfg.std_file}")
    info(f"wght_out={cfg.bio_weight_file}")
    info(f"ogcm(std)={ogcm_file}")
    info(f"bio_model={cfg.bio_model_type}")

# [02] Load bio_vars.yml (ini vars + rules)
with step("[02] Load bio vars (YAML)"):
    bio_yaml = getattr(cfg, "bio_yaml", "bio_vars.yml")
    bio_db = load_bio_yaml(bio_yaml)

    bio_model = str(cfg.bio_model_type).lower()

    # NOTE: std uses the same variable list and rules as ini
    ini_vars = get_bio_vars(bio_db, bio_model, target="ini")
    bio_rules = get_bio_rules(bio_db, bio_model, target="ini")

    info(f"bio_yaml={bio_yaml}")
    info(f"ini_vars={ini_vars}")

# [03] Crop OGCM domain and prepare remap weights
with step("[03] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.bio_weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.bio_weight_file}")
        plus(f"Weight file created: {cfg.bio_weight_file}")
    else:
        info(f"Use existing wght file {cfg.bio_weight_file}")

# [04] Load weights once
with step("[04] Load remap weights"):
    with Dataset(cfg.bio_weight_file) as ncw:
        row = ncw.variables["row"][:] - 1
        col = ncw.variables["col"][:] - 1
        S   = ncw.variables["S"][:]

note(f"Flood method: {cfg.flood_method_for_ini}")

# [05] Open OGCM std file once
nc_raw = Dataset(ogcm_file, mode="r", maskandscale=True)
nc_clim = tl.MaskedNetCDF(nc_raw)

try:
    for idt in range(12):

        mtag = month_tag(idt, cfg)
        outname = os.path.join(cfg.std_dir, f"{cfg.std_file}_{mtag}.nc")

        if not os.path.exists(outname):
            raise FileNotFoundError(f"STD output file not found (create_std first): {outname}")

        # [05] Load OGCM fields (yaml-driven) + rules
        with step("[05] Load OGCM fields", ts=f"idt={idt}"):
            field_dict = {}

            # (1) OGCM에서 직접 읽을 수 있는 것만 로딩
            for var in ini_vars:
                src = bio_map.get(var, None)
                if src is None:
                    continue
                field_dict[var] = nc_clim.get(src, idt, slice(None), idy, idx)


            # >>> FIXME >>>
            print("!!! SHJO: FIX ME !!!")
            factor = (0.02 * 6.625 * 12.0)
            val = field_dict["phytoplankton"]

            if isinstance(val, np.ma.MaskedArray):
                val = val.filled(np.nan)  # 마스크 → NaN 유지
            field_dict["phytoplankton"] = val / factor

#            field_dict["phytoplankton"]=field_dict["phytoplankton"]/(0.02*6.625*12) 
#            val = field_dict["phytoplankton"]
#            print("BEFORE scale", np.nanmin(val), np.nanmax(val), np.isnan(val).sum())
#            val2 = val / (0.02*6.625*12)
#            print("AFTER  scale", np.nanmin(val2), np.nanmax(val2), np.isnan(val2).sum())
#            field_dict["phytoplankton"] = val2
#            print("!!!")
#            val = field_dict["phytoplankton"]
#            print(type(val), isinstance(val, np.ma.MaskedArray))
#            print("fill_value:", getattr(val, "fill_value", None))
#            v = np.array(val)  # 강제로 ndarray
#            print("ndarray nan:", np.isnan(v).sum(), "minmax:", np.nanmin(v), np.nanmax(v))
#            print("!!!")
            # <<< FIXME <<<

            # (2) rules로 파생/상수 채우기
            apply_bio_rules(field_dict, bio_rules)

            # (3) 누락 체크
            missing = [v for v in ini_vars if v not in field_dict]
            if missing:
                raise KeyError(f"bio missing vars (no mapping & no rule): {missing}")

            field = tl.ConfigObject(**field_dict)

        # [06] Apply remap weights to all fields
        with step("[06] Remap (weights)"):
            for varname in ini_vars:
                var_src = getattr(field, varname)
                with capture_warnings(tag="remap"):
                    remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
                setattr(field, varname, remapped)

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
        with step("[09] Mask & clean"):
            for var in ini_vars:
                arr = getattr(field, var)
                if isinstance(arr, np.ndarray) and arr.ndim == 2:
                    arr[grd.mask_rho == 0] = 0.0
                else:
                    arr[..., grd.mask_rho == 0] = 0.0
                setattr(field, var, arr)

        # [10] z→σ interpolation (bio: rho-point only)
        with step("[10] z→σ interpolation"):
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

            for var in ini_vars:
                val = getattr(field, var)

                # 3D only
                if not isinstance(val, np.ndarray) or val.ndim != 3:
                    continue

                padded = np.vstack((val[0:1], val, val[-1:]))
                flipped = np.flip(padded, axis=0)
                with capture_warnings(tag="z2sigma"):
                    remapped = tl.ztosigma_numba(flipped, zr, Zf)
                setattr(field, var, remapped)

        # [11] Write variables to std file
        with step("[11] Write variables", out=ellipsis(outname)):
            with Dataset(outname, 'a') as nc:
                for varname in ini_vars:
                    if varname not in nc.variables:
                        warn_line(f"skip write (not in std): {varname}")
                        continue
                    nc[varname][0] = getattr(field, varname) * cfg.bio_scale_factor

finally:
    nc_raw.close()

bar("Summary")
print(f"Total elapsed: {time.time() - start0:.3f}s")
