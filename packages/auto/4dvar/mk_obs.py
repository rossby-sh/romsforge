
# mk_obs.py
# - mk_bry 스타일로 time 수집: collect_time_info(filelist, time_var, (start,end))
# - survey dimension: fixed (A) = len(tinfo)
# - datum dimension: unlimited (append)
# - config key 통일: longitude/latitude 사용 (파일이 lon/lat이면 config에서 매핑)
# - cfg.is3d == False -> obs_depth = 0.0

# --- [00] Imports and path setup ---
import sys
import os
import glob
import time
from collections import defaultdict

import numpy as np
from netCDF4 import Dataset

# libs path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..", "libs")))
import create as cn
import utils as tl
from io_utils import collect_time_info  # <-- legacy 말고 range 버전

# logging
from log_utils2 import configure, step, capture_warnings, info, warn_line, bar
configure(width=80, show_sections=False, color_mode="auto")


# --- [helpers] ---------------------------------------------------------------
def _read_time_units(nc_path: str, time_name: str) -> str:
    with Dataset(nc_path) as nc:
        tvar = nc.variables[time_name]
        tunits = getattr(tvar, "units", None)
        if tunits is None:
            raise RuntimeError(f"time variable has no units attribute: {nc_path} ({time_name})")
    return tunits


def build_obs_chunk_from_source(nc_path, idt, rel_time, grd, cfg):
    """
    Build 1D datum-chunk arrays from a single source file/time index.

    Required cfg:
      cfg.src_var_name: {"time":..., "chlor_a":..., "longitude":..., "latitude":...}
      cfg.thin_y, cfg.thin_x
      cfg.valid_min, cfg.valid_max
      cfg.unit_conv
      cfg.err_rel, cfg.err_floor
      cfg.obs_flag, cfg.obs_provenance
      cfg.is3d (bool), cfg.obs_depth_const (used only when is3d=True)
    """
    with Dataset(nc_path, maskandscale=True) as nc_raw:
        nc = tl.MaskedNetCDF(nc_raw)

        # chlor_a: expect (time, lat, lon) or (lat, lon)
        chl = nc.get(cfg.src_var_name["chlor_a"], idt, slice(None), slice(None))

        # lon/lat: MODIS는 보통 1D (lon(lon), lat(lat))
        lon = nc.get(cfg.src_var_name["longitude"], slice(None))
        lat = nc.get(cfg.src_var_name["latitude"],  slice(None))

    chl = np.asarray(chl)
    lon = np.asarray(lon)
    lat = np.asarray(lat)

    # lon/lat grid -> 2D
    if lon.ndim == 1 and lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    elif lon.ndim == 2 and lat.ndim == 2:
        lon2d, lat2d = lon, lat
    else:
        raise RuntimeError(
            f"Unexpected lon/lat dims in {os.path.basename(nc_path)}: "
            f"lon.ndim={lon.ndim}, lat.ndim={lat.ndim}"
        )

    # thinning
    ty, tx = int(cfg.thin_y), int(cfg.thin_x)
    chl   = chl[::ty, ::tx]
    lon2d = lon2d[::ty, ::tx]
    lat2d = lat2d[::ty, ::tx]

    # valid filter
    valid = np.isfinite(chl) & (chl > cfg.valid_min) & (chl < cfg.valid_max)
    if not np.any(valid):
        return None

    obs_val = chl[valid] * float(cfg.unit_conv)
    obs_lon = lon2d[valid]
    obs_lat = lat2d[valid]

    Nobs = int(obs_val.size)
    if Nobs == 0:
        return None

    # time/type/prov
    obs_time = np.full(Nobs, float(rel_time), dtype=np.float64)
    obs_type = np.full(Nobs, int(cfg.obs_flag), dtype=np.int32)
    obs_prov = np.full(Nobs, float(cfg.obs_provenance), dtype=np.float64)

    # depth policy (너 확정)
    if bool(getattr(cfg, "is3d", False)):
        depth_val = float(cfg.obs_depth_const)
    else:
        depth_val = 0.0
    obs_dep = np.full(Nobs, depth_val, dtype=np.float64)
    obs_zg  = np.zeros(Nobs, dtype=np.float64)

    # error policy
    obs_err = np.maximum((obs_val * float(cfg.err_rel)) ** 2, float(cfg.err_floor)).astype(np.float64)

    # grid position: use obs_ijpos
#    with capture_warnings(tag="obs_ijpos"):
#        Xg, Yg = tl.obs_ijpos(obs_lon, obs_lat, grd)
#    obs_Xg = np.asarray(Xg, dtype=np.float64)
#    obs_Yg = np.asarray(Yg, dtype=np.float64)
    obs_Xg,obs_Yg=obs_lon,obs_lat

    return {
        "Nobs": Nobs,
        "obs_type": obs_type,
        "obs_time": obs_time,
        "obs_lon": obs_lon.astype(np.float64),
        "obs_lat": obs_lat.astype(np.float64),
        "obs_depth": obs_dep,
        "obs_error": obs_err,
        "obs_value": obs_val.astype(np.float64),
        "obs_Xgrid": obs_Xg,
        "obs_Ygrid": obs_Yg,
        "obs_Zgrid": obs_zg,
        "obs_provenance": obs_prov,
    }


# --- main --------------------------------------------------------------------
start0 = time.time()
bar("Observation Build (obs)")

# [01] Load configuration and input metadata
with step("[01] Load configuration and input metadata"):
    cfg = tl.parse_config("./config_obs.yaml")
    grd = tl.load_roms_grid(cfg.grdname)

    info(f"grid={cfg.grdname}")
    info(f"obs_out={cfg.obs_out}")
    info(f"src_glob={cfg.obs_src_glob}")

# [02] Time index matching & relative time (mk_bry 스타일)
with step("[02] Time index matching & relative time"):
    filelist = sorted(glob.glob(cfg.obs_src_glob))
    if not filelist:
        raise RuntimeError(f"No files matched: {cfg.obs_src_glob}")

    tinfo = collect_time_info(
        filelist,
        cfg.src_var_name["time"],
        (str(cfg.obs_start_date), str(cfg.obs_end_date)),
    )
    if not tinfo:
        raise RuntimeError("No valid time steps found in given range.")

    datenums = np.array([ti.raw_value for ti in tinfo], dtype=np.float64)

    time_unit = _read_time_units(tinfo[0].filename, cfg.src_var_name["time"])
    relative_time = tl.compute_relative_time(datenums, time_unit, cfg.time_ref).astype(np.float64)

    # survey_time은 mk_bry처럼 relative_time을 그대로 사용
    survey_time = np.asarray(relative_time, dtype=np.float64)

    # seq 붙여서 grouping해도 순서 유지
    tinfo_seq = list(enumerate(tinfo))  # (seq, TimeIndex)

    info(f"Ns = {len(tinfo_seq)}")
    info(f"time_units={time_unit}")

# [03] List & group source files (open 최소화)
with step("[03] List & group source files"):
    grouped = defaultdict(list)  # filename -> list[(seq, TimeIndex)]
    for seq, entry in tinfo_seq:
        grouped[entry.filename].append((seq, entry))
    for entries in grouped.values():
        entries.sort(key=lambda x: x[1].datetime)

# [04] Create obs NetCDF skeleton (datum=unlimited, survey=fixed)
with step("[04] Create obs NetCDF skeleton", out=cfg.obs_out):
    if getattr(cfg, "force_write", False) and os.path.exists(cfg.obs_out):
        os.remove(cfg.obs_out)

    status = cn.create_obs(
        cfg,
        grd,
        outname=cfg.obs_out,
        survey_time=survey_time,
        Nsurvey=len(tinfo_seq),        # A: fixed
        ncFormat=cfg.ncformat,
        unlimited_datum=True,          # datum=None (append)
    )
    if status:
        raise RuntimeError(f"Failed creating obs file: {cfg.obs_out}")

# [05] File loop -> build chunk -> append (survey별 Nobs 기록)
with step("[05] Build & append obs chunks"):
    Nobs_by_survey = np.zeros(len(tinfo_seq), dtype=np.int32)

    for filename, entries in grouped.items():
        with step("[05] Open source file", file=os.path.basename(filename)):
            for seq, e in entries:
                rel_time = float(survey_time[seq])
                idt = int(e.index)

                chunk = build_obs_chunk_from_source(filename, idt, rel_time, grd, cfg)
                if chunk is None:
                    Nobs_by_survey[seq] = 0
                    continue

                Nobs_by_survey[seq] = int(chunk["Nobs"])
                cn.append_obs(cfg.obs_out, chunk, ncFormat=cfg.ncformat)

# [06] Finalize survey vars + variance
with step("[06] Finalize survey table"):
    cn.finalize_obs(
        cfg.obs_out,
        survey_time=survey_time,
        Nobs=Nobs_by_survey,
        obs_variance=np.asarray(cfg.obs_variance, dtype=np.float64),
        ncFormat=cfg.ncformat,
    )

bar("Summary")
print(f"Total elapsed: {time.time() - start0:.3f}s")
