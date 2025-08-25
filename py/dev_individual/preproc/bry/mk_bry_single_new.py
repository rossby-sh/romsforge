# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
from collections import defaultdict
import time
import contextlib
import warnings
import traceback

# Append libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_B as cn
import utils as tl
from io_utils import collect_time_info

# === pretty log helpers =======================================================
TERMW = 80  # adjust to 80~100 by preference

def bar(title: str):
    t = f"== {title} "
    fill = "=" * max(0, TERMW - len(t))
    print(t + fill)

def sec(title: str):
    t = f"-- {title} "
    fill = "-" * max(0, TERMW - len(t))
    print(t + fill)

def ellipsis(path: str, keep: int = 2) -> str:
    # /a/b/c/d/e.nc -> /.../d/e.nc
    parts = path.split(os.sep)
    if len(parts) <= keep + 1:
        return path
    return os.sep + "..." + os.sep + os.sep.join(parts[-keep:])

def done(ts: str, dur: float):
    # [DONE] 2025-05-01 00 | 3.454s
    print(f"[DONE] {ts:<13} | {dur:>6.3f}s")

def info(line: str):
    print(f"· {line}")

def warn_line(msg: str):
    print(f"[WARN] {msg}")

@contextlib.contextmanager
def step(title: str, **kv):
    """Pretty step block. On exception, prints FAIL with brief location and re-raises."""
    meta = " ".join(f"{k}={v}" for k, v in kv.items()) if kv else ""
    head = (title if not meta else f"{title} | {meta}")
    #sec(head)
    t0 = time.time()
    try:
        yield
    except Exception as e:
        el = time.time() - t0
        print(f"[FAIL] {title} | dur={el:.3f}s | {type(e).__name__}: {e}")
        tb = traceback.extract_tb(e.__traceback__)
        if tb:
            src = tb[-1]
            print(f"       at {os.path.basename(src.filename)}:{src.lineno} in {src.name}")
        raise
    else:
        el = time.time() - t0
        print(f"[OK]   {title} | dur={el:.3f}s")

@contextlib.contextmanager
def capture_warnings(tag: str | None = None):
    """Capture numpy/python warnings and print them as pretty [WARN] lines after the block."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        yield
        for wi in w:
            fn = os.path.basename(getattr(wi, "filename", "?"))
            ln = getattr(wi, "lineno", "?")
            msg = str(wi.message)
            lbl = f"{fn}:{ln}"
            warn_line(f"{lbl} — {msg}" + (f" | {tag}" if tag else ""))

# numpy warnings: show as warnings (not silent)
np.seterr(all="warn")
# ============================================================================

# --- main --------------------------------------------------------------------
start1 = time.time()
bar("Boundary Build")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg = tl.parse_config("./config_single.yaml")
    grd = tl.load_roms_grid(cfg.grdname)
    filelist = sorted(glob.glob(os.path.join(cfg.ogcm_path, "*.nc")))
    assert len(filelist) > 0, "No OGCM files found in ogcm_path!"
    ogcm = tl.load_ogcm_metadata(filelist[0], cfg.ogcm_var_name)
    info(f"src=OGCM grid={os.path.basename(cfg.grdname)}")
    info(f"out={ellipsis(cfg.bryname)} wght={ellipsis(cfg.weight_file)}")

# [02] Time index matching and relative time calculation
with step("[02] Time index matching & relative time"):
    tinfo = collect_time_info(
        filelist,
        cfg.ogcm_var_name.time,
        (str(cfg.bry_start_date), str(cfg.bry_end_date)),
    )
    datenums = np.array([ti.raw_value for ti in tinfo])
    relative_time = tl.compute_relative_time(datenums, ogcm.time_unit, cfg.time_ref)

# [03] Create initial NetCDF file
with step("[03] Create initial NetCDF"):
    status = cn.create_bry(
        cfg, grd, relative_time, bio_model=cfg.bio_model_type, ncFormat=cfg.ncformat
    )
    if status:
        raise RuntimeError("create_bry() returned error")

# [04] Crop OGCM domain and prepare remap weights
with step("[04] Prepare weights", reuse=not cfg.calc_weight):
    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(
        ogcm.lat, ogcm.lon, grd.lat, grd.lon
    )
    if cfg.calc_weight:
        status = tl.build_bilinear_regridder(
            lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False
        )
        if status:
            raise RuntimeError("build_bilinear_regridder() failed")
    else:
        info(f"Use existing wght file {ellipsis(cfg.weight_file)}")
    with Dataset(cfg.weight_file) as nc:
        row = nc.variables["row"][:] - 1
        col = nc.variables["col"][:] - 1
        S = nc.variables["S"][:]

sec("Notes")
info(f"Biological variables: {cfg.bio_model_type}")
info(f"Flood method (boundary): {cfg.flood_method_for_bry}")

# Prep bry data containers
bry_data = tl.make_all_bry_data_shapes(
    ["zeta", "ubar", "vbar", "temp", "salt", "u", "v"],
    len(tinfo),
    grd,
    cfg.vertical.layer_n,
)

bry_time = []

# [05] List & group OGCM files
with step("[05] List & group OGCM files"):
    grouped = defaultdict(list)
    for entry in tinfo:
        grouped[entry.filename].append(entry)
    for entries in grouped.values():
        entries.sort(key=lambda x: x.datetime)
    time_index_map = {entry.datetime: n for n, entry in enumerate(tinfo)}

# --- File loop ---------------------------------------------------------------
for filename, entries in grouped.items():
    with step("[06] Open source file", file=os.path.basename(filename)):
        start2 = time.time()
        with Dataset(filename, maskandscale=True) as nc:
            nc_wrap = tl.MaskedNetCDF(nc)

            for entry in entries:
                i = entry.index
                t = entry.datetime
                tval = entry.raw_value
                n = time_index_map[t]
                tag = f"ts={t:%Y-%m-%d %H} file={os.path.basename(filename)}"

                # [07] Load OGCM raw fields
                with step("[07] Load OGCM fields", ts=f"{t:%Y-%m-%d %H}"):
                    zeta = nc_wrap.get(cfg.ogcm_var_name['zeta'], i, idy, idx)
                    temp = nc_wrap.get(cfg.ogcm_var_name['temperature'], i, slice(None), idy, idx)
                    salt = nc_wrap.get(cfg.ogcm_var_name['salinity'], i, slice(None), idy, idx)
                    u    = nc_wrap.get(cfg.ogcm_var_name['u'], i, slice(None), idy, idx)
                    v    = nc_wrap.get(cfg.ogcm_var_name['v'], i, slice(None), idy, idx)
                    with capture_warnings(tag):
                        ubar = tl.depth_average(u, ogcm.depth)
                        vbar = tl.depth_average(v, ogcm.depth)
                    field = tl.ConfigObject(zeta=zeta, ubar=ubar, vbar=vbar,
                                            temp=temp, salt=salt, u=u, v=v)

                # [08] Remap (weights)
                with step("[08] Remap (weights)", ts=f"{t:%Y-%m-%d %H}"):
                    for var in vars(field):
                        var_src = getattr(field, var)
                        with capture_warnings(tag):
                            remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
                        setattr(field, var, remapped)

                # [09] Flood H/V
                with step("[09] Flood H/V", ts=f"{t:%Y-%m-%d %H}"):
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

                # [10] Mask & rotate
                with step("[10] Mask & rotate", ts=f"{t:%Y-%m-%d %H}"):
                    for var in vars(field):
                        arr = getattr(field, var)
                        if arr.ndim == 2:
                            arr[grd.mask == 0] = 0.0
                        else:
                            arr[..., grd.mask == 0] = 0.0
                        setattr(field, var, arr)
                    u_rot, v_rot       = tl.rotate_vector_euler(field.u,    field.v,    grd.angle, to_geo=False)
                    ubar_rot, vbar_rot = tl.rotate_vector_euler(field.ubar, field.vbar, grd.angle, to_geo=False)
                    setattr(field, 'u',    tl.rho2uv(u_rot, 'u'))
                    setattr(field, 'v',    tl.rho2uv(v_rot, 'v'))
                    setattr(field, 'ubar', tl.rho2uv(ubar_rot, 'u'))
                    setattr(field, 'vbar', tl.rho2uv(vbar_rot, 'v'))

                # [11] z→σ & save bry
                with step("[11] z→σ & save bry", ts=f"{t:%Y-%m-%d %H}"):
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

                    zr_3d = tl.zlevs(*zargs, 1, grd.topo, field.zeta)    # (Ns, Mp, Lp)
                    zw_3d = tl.zlevs(*zargs, 5, grd.topo, field.zeta)    # (Ns+1, Mp, Lp)
                    dz_3d = zw_3d[1:, :, :] - zw_3d[:-1, :, :]

                    zu_3d = tl.rho2uv(zr_3d, "u")
                    zv_3d = tl.rho2uv(zr_3d, "v")
                    dz_u3d = tl.rho2uv(dz_3d, "u")
                    dz_v3d = tl.rho2uv(dz_3d, "v")

                    for direction in directions:
                        zr = tl.extract_bry(zr_3d, direction)
                        zu = tl.extract_bry(zu_3d, direction)
                        zv = tl.extract_bry(zv_3d, direction)
                        dzr = tl.extract_bry(dz_3d, direction)
                        dzu = tl.extract_bry(dz_u3d, direction)
                        dzv = tl.extract_bry(dz_v3d, direction)

                        for var, zgrid in zip(['temp', 'salt', 'u', 'v'], [zr, zr, zu, zv]):
                            val = tl.extract_bry(getattr(field, var), direction)        # (Nz, Lp)
                            val_pad = np.vstack((val[0:1], val, val[-1:]))              # (Nz+2, Lp)
                            val_flip = np.flip(val_pad, axis=0)                         # (Nz+2)
                            with capture_warnings(tag):
                                sigma_interped = tl.ztosigma_1d_numba(val_flip, zgrid, Z_flipped)  # (Ns, Lp)

                            if var == 'u':
                                barotropic = tl.extract_bry(field.ubar, direction)
                                dz_bar = dzu
                            elif var == 'v':
                                barotropic = tl.extract_bry(field.vbar, direction)
                                dz_bar = dzv
                            else:
                                barotropic = None

                            if barotropic is not None:
                                with capture_warnings(tag):
                                    sigma_interped, _, barotropic_corrected, _ = tl.conserve_and_recompute_barotropic(
                                        sigma_interped, sigma_interped, barotropic, barotropic, dz_bar, dz_bar
                                    )
                                setattr(field, f"{'ubar' if var == 'u' else 'vbar'}_{direction}", barotropic_corrected)

                            setattr(field, f"{var}_{direction}", sigma_interped)

                        # zeta separate
                        val = tl.extract_bry(field.zeta, direction)
                        setattr(field, f"zeta_{direction}", val)

                        # save to buffers
                        for varname in bry_data:
                            val_d = getattr(field, f"{varname}_{direction}")
                            bry_data[varname][direction][n, ...] = val_d

                        # time save
                        time_converted = tl.compute_relative_time(tval, ogcm.time_unit, cfg.time_ref)
                        bry_time.append(time_converted)

                # per-timestep done line
                done(f"{t:%Y-%m-%d %H}", time.time() - start2)

# [12] Write all remapped variables to bry.nc
with step("[12] Write variables", out=ellipsis(cfg.bryname)):
    with Dataset(cfg.bryname, 'a') as nc:
        # If time variables are needed later, uncomment and map accordingly
        # for tname in ['bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']:
        #     nc[tname][:] = np.array(bry_time)
        for varname in bry_data:
            for direction in bry_data[varname]:
                var_fullname = f"{varname}_{direction}"
                nc[var_fullname][:] = bry_data[varname][direction]

bar("Summary")
print(f"Total elapsed: {time.time() - start1:.3f}s")

