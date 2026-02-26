
# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
from netCDF4 import Dataset
from collections import defaultdict
import time
import datetime as dt

# Append libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl
from io_utils import collect_time_info, load_bio_yaml, get_bio_vars, get_bio_rules, apply_bio_rules
# logging
from log_utils2 import configure, step, capture_warnings, info, note, plus, bar, warn_line, done, ellipsis
configure(width=80, show_sections=False, color_mode='auto')

# --- main --------------------------------------------------------------------
start1 = time.time()
bar("Boundary Build (bio)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg = tl.parse_config("./config.yaml")
    grd = tl.load_roms_grid(cfg.grdname)

    ogcm_dir = cfg.bio_ogcm_name
    bio_map = dict(cfg.bio_ogcm_var_name.__dict__)

    assert os.path.isdir(ogcm_dir), f"OGCM daily directory not found: {ogcm_dir}"

    filelist = sorted(glob.glob(os.path.join(ogcm_dir, "cmems_bio_*.nc")))
    assert len(filelist) > 0, f"No daily CMEMS files found in: {ogcm_dir}"

    info(f"src=OGCM(daily dir) grid={os.path.basename(cfg.grdname)}")
    info(f"out={ellipsis(cfg.bryname)} wght={ellipsis(cfg.bio_weight_file)}")
    info(f"ogcm_dir={ogcm_dir}")
    info(f"daily_files={len(filelist)}")

# [02] Load bio_vars.yml (bry vars + rules)
with step("[02] Load bio vars (YAML)"):
    bio_yaml = getattr(cfg, "bio_yaml", "bio_vars.yml")
    bio_db = load_bio_yaml(bio_yaml)

    bry_vars = get_bio_vars(bio_db, str(cfg.bio_model_type).lower(), target="bry")

    bio_rules = get_bio_rules(bio_db, cfg.bio_model_type,target="bry")

    info(f"bio_yaml={bio_yaml}")
    info(f"bio_model={cfg.bio_model_type}")
    info(f"bry_vars={bry_vars}")

# [03] Time index matching & relative time
with step("[03] Time index matching & relative time"):
    def _ymd(x):
        return x[:10] if isinstance(x, str) else f"{x:%Y-%m-%d}"

    #time_var = cfg.bio_ogcm_var_name["time"]

    time_var = bio_map["time"]
    tinfo = collect_time_info(
        filelist,
        time_var,
        (_ymd(cfg.bry_start_date), _ymd(cfg.bry_end_date))
    )
    info(f"bry steps={len(tinfo)}")

# [04] Prepare weights (once)
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    ogcm0 = tl.load_ogcm_metadata(tinfo[0].filename, cfg.bio_ogcm_var_name)
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm0.lat, ogcm0.lon, grd.lat, grd.lon)

    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.bio_weight_file, reuse=False)
        if status:
            raise RuntimeError(f"Failed to generate remap weights: {cfg.bio_weight_file}")
        plus(f"Weight file created: {cfg.bio_weight_file}")
    else:
        info(f"Use existing wght file {cfg.bio_weight_file}")

    with Dataset(cfg.bio_weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S   = nc.variables["S"][:]

# [05] Allocate bry buffers
with step("[05] Allocate bry buffers"):
    bry_data = tl.make_all_bry_data_shapes(
        bry_vars,
        len(tinfo),
        grd,
        cfg.vertical.layer_n,
    )

# [06] Group time entries by file
with step("[06] Group time entries by file"):
    grouped = defaultdict(list)
    for entry in tinfo:
        grouped[entry.filename].append(entry)
    for entries in grouped.values():
        entries.sort(key=lambda x: x.datetime)
    time_index_map = {entry.datetime: n for n, entry in enumerate(tinfo)}

note(f"Flood method: {cfg.flood_method_for_bry}")

# [07] Main loop
for filename, entries in grouped.items():
    with Dataset(filename, maskandscale=True) as nc:
        nc_wrap = tl.MaskedNetCDF(nc)
        ogcm = tl.load_ogcm_metadata(filename, cfg.bio_ogcm_var_name)

        for entry in entries:
            i = entry.index
            t = entry.datetime
            n = time_index_map[t]
            tag = f"{t:%Y-%m-%d}"

            # [07-1] Load OGCM fields (mapping-driven) + rules
            with step("[07] Load OGCM fields", ts=tag):
                field_dict = {}

                for var in bry_vars:
                    if var in bio_map:
                        src = bio_map[var]
                        field_dict[var] = nc_wrap.get(src, i, slice(None), idy, idx)

                # >>> FIXME >>>
                print("!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!")
                field_dict["phyt"]=field_dict["phyt"]/(0.02*6.625*12)
                # <<< FIXME <<<
                apply_bio_rules(field_dict, bio_rules)

                missing = [v for v in bry_vars if v not in field_dict]
                if missing:
                    raise KeyError(f"bry missing vars (no mapping & no rule): {missing}")

                field = tl.ConfigObject(**field_dict)

            # [07-2] Remap
            with step("[08] Remap (weights)", ts=tag):
                for var in bry_vars:
                    var_src = getattr(field, var)
                    with capture_warnings(tag="remap"):
                        remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
                    setattr(field, var, remapped)

            # [07-3] Flood: horizontal
            with step("[09] Flood: horizontal", ts=tag):
                for var in bry_vars:
                    val = getattr(field, var)
                    val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_bry)
                    setattr(field, var, val_flooded)

            # [07-4] Flood: vertical
            with step("[10] Flood: vertical", ts=tag):
                for var in bry_vars:
                    val = getattr(field, var)
                    if not isinstance(val, np.ndarray) or val.ndim != 3:
                        continue
                    with capture_warnings(tag="vflood"):
                        val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask_rho), spval=-1e10)
                    setattr(field, var, val_flooded)

            # [07-5] Mask
            with step("[11] Mask", ts=tag):
                for var in bry_vars:
                    arr = getattr(field, var)
                    if arr.ndim == 2:
                        arr[grd.mask_rho == 0] = 0.0
                    else:
                        arr[..., grd.mask_rho == 0] = 0.0
                    setattr(field, var, arr)

            # [07-6] z→σ & save bry (bio vars are all rho-point → zr)
            with step("[12] z→σ & save", ts=tag):
                directions = ["north", "south", "east", "west"]

                zargs = (
                    cfg.vertical.vtransform,
                    cfg.vertical.vstretching,
                    cfg.vertical.theta_s,
                    cfg.vertical.theta_b,
                    cfg.vertical.tcline,
                    cfg.vertical.layer_n,
                )
                zr_3d = tl.zlevs(*zargs, 1, grd.topo, np.zeros_like(grd.topo))

                Z = np.zeros(len(ogcm.depth) + 2)
                Z[0] = 100
                Z[1:-1] = -np.abs(ogcm.depth)
                Z[-1] = -100000
                Z_flipped = np.flipud(Z)

                for direction in directions:
                    zr = tl.extract_bry(zr_3d, direction)

                    for var in bry_vars:
                        val = tl.extract_bry(getattr(field, var), direction)

                        # 2D (surface)면 z→σ 없이 그대로 쓴다
                        if not isinstance(val, np.ndarray) or val.ndim == 1:
                            bry_data[var][direction][n, ...] = val
                            continue

                        # 3D (z, Lp) → sigma
                        val_pad = np.vstack((val[0:1], val, val[-1:]))
                        val_flip = np.flip(val_pad, axis=0)

                        with capture_warnings(tag="z2sigma"):
                            sigma = tl.ztosigma_1d_numba(val_flip, zr, Z_flipped)

                        bry_data[var][direction][n, ...] = sigma

# [13] Write variables to bry.nc
with step("[13] Write variables", out=ellipsis(cfg.bryname)):
    with Dataset(cfg.bryname, "a") as nc:
        for var in bry_vars:
            for direction in ("north", "south", "east", "west"):
                vname = f"{var}_{direction}"
                if vname not in nc.variables:
                    warn_line(f"skip write (not in bry): {vname}")
                    continue
                nc[vname][:] = bry_data[var][direction]

bar("Summary")
print(f"--- Time elapsed: {time.time()-start1:.3f}s ---")
