

# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import create_B as cn
import utils as tl
from io_utils import collect_time_info
from collections import defaultdict
import time 

from numba import set_num_threads
set_num_threads(4)

# --- [01] Load configuration and input metadata ---
cfg  = tl.parse_config("./config_single.yaml")
grd  = tl.load_roms_grid(cfg.grdname)
filelist = sorted(glob.glob(os.path.join(cfg.ogcm_path, "*.nc")))

assert len(filelist) > 0, "No HYCOM files found in ogcm_path!"

ogcm = tl.load_ogcm_metadata(filelist[0], cfg.ogcm_var_name)

start1=time.time()
# --- [02] Time index matching and relative time calculation ---
tinfo = collect_time_info(filelist, cfg.ogcm_var_name['time'],\
        (str(cfg.bry_start_date), str(cfg.bry_end_date)) )
print(f'--- Time elapsed: {time.time()-start1:.3f}s ---')

datenums = np.array([tval for _, _, _, tval in tinfo])

relative_time = tl.compute_relative_time(datenums, ogcm.time_unit, cfg.time_ref)

# --- [03] Create initial NetCDF file ---
print(f"--- [03] Creating boundary NetCDF file: {cfg.bryname} ---")
status = cn.create_bry(cfg, grd, relative_time,  bio_model=cfg.bio_model_type, ncFormat=cfg.ncformat)
if status:
    print(f"--- [!ERROR] Failed to creating file {cfg.bryname} ---")
    raise
print(f"--- Created file: {cfg.bryname} ---")

# --- [04] Crop OGCM domain and prepare remap weights ---
lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)

if cfg.calc_weight:
    print(f"--- [04] Calculating weight: {cfg.weight_file} ---")
    status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
    
    if status:
        print(f"--- [!ERROR] Failed to generate remap weights: {cfg.weight_file} ---")
        raise  
else: 
    print(f"--- [04] Use existing wght file {cfg.weight_file} ---")



with Dataset(cfg.weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]

bry_data = {
    'zeta': {d: [] for d in ['west', 'east', 'south', 'north']},
    'temp': {d: [] for d in ['west', 'east', 'south', 'north']},
    'salt': {d: [] for d in ['west', 'east', 'south', 'north']},
    'u':    {d: [] for d in ['west', 'east', 'south', 'north']},
    'v':    {d: [] for d in ['west', 'east', 'south', 'north']},
    'ubar': {d: [] for d in ['west', 'east', 'south', 'north']},
    'vbar': {d: [] for d in ['west', 'east', 'south', 'north']}
}

bry_time = []

# --- [05] Load OGCM raw fields (zeta, temp, salt, u, v) ---
print("--- [05] Loading OGCM data ---")
grouped = defaultdict(list)
for f, i, t, tval in tinfo:
    grouped[f].append((i, t, tval))

# --- 정렬 추가 ---
for f in grouped:
    grouped[f].sort(key=lambda x: x[1])  # datetime 기준 정렬

import time

# --- 파일 단위 루프 ---
for f, entries in grouped.items():
    with Dataset(f, maskandscale=True) as nc:
        nc_wrap = tl.MaskedNetCDF(nc)
        start=time.time()
        for i, t, tval in entries:
            # [05] OGCM 필드 로딩
            zeta = nc_wrap.get(cfg.ogcm_var_name['zeta'], i, idy, idx)
            temp = nc_wrap.get(cfg.ogcm_var_name['temperature'], i, slice(None), idy, idx)
            salt = nc_wrap.get(cfg.ogcm_var_name['salinity'], i, slice(None), idy, idx)
            u    = nc_wrap.get(cfg.ogcm_var_name['u'], i, slice(None), idy, idx)
            v    = nc_wrap.get(cfg.ogcm_var_name['v'], i, slice(None), idy, idx)
            ubar = tl.depth_average(u, ogcm.depth)
            vbar = tl.depth_average(v, ogcm.depth)

            field = tl.ConfigObject(zeta=zeta, ubar=ubar, vbar=vbar,
                                    temp=temp, salt=salt, u=u, v=v)

            # [06] Remap
            print("--- [06] Remapping ---")
            for var in vars(field):
                var_src = getattr(field, var)
                remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
                setattr(field, var, remapped)

            # [07] Horizontal Flood
            print("--- [07] Apply horizontal flood ---")
            for var in vars(field):
                val = getattr(field, var)
                val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_bry)
                setattr(field, var, val_flooded)

            print("--- [08] No vertical flood ---")
            for var in vars(field):
                val = getattr(field, var)
                if var.ndim == 2:
                    continue
                #val_flooded = tl.flood_vertical_vectorized(val, grd.mask, spval=-1e10)
                val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask), spval=-1e10)
                setattr(field, var, val_flooded)

            # [09] Mask land
            print("--- [09] Masking land to 0 ---")
            for var in vars(field):
                arr = getattr(field, var)
                if arr.ndim==2:
                    arr[grd.mask == 0] = 0.0
                else:
                    arr[...,grd.mask == 0] = 0.0
                setattr(field, var, arr)


            # [10] Rotate
            print("--- [10] Rotate vectors ---")
            u_rot, v_rot       = tl.rotate_vector_euler(field.u,    field.v,    grd.angle, to_geo=False)
            ubar_rot, vbar_rot = tl.rotate_vector_euler(field.ubar, field.vbar, grd.angle, to_geo=False)
            setattr(field, 'u',    tl.rho2uv(u_rot, 'u'))
            setattr(field, 'v',    tl.rho2uv(v_rot, 'v'))
            setattr(field, 'ubar', tl.rho2uv(ubar_rot, 'u'))
            setattr(field, 'vbar', tl.rho2uv(vbar_rot, 'v'))

            # [11] Vertical interpolation
            print("--- [11] Vertical interpolation from z to sigma ---")
            zargs = (
                cfg.vertical.vtransform,
                cfg.vertical.vstretching,
                cfg.vertical.theta_s,
                cfg.vertical.theta_b,
                cfg.vertical.tcline,
                cfg.vertical.layer_n,
            )

            zr = tl.zlevs(*zargs, 1, grd.topo, field.zeta)
            zu, zv = tl.rho2uv(zr, 'u'), tl.rho2uv(zr, 'v')

            Z = np.zeros(len(ogcm.depth) + 2)
            Z[0] = 100; Z[1:-1] = -np.abs(ogcm.depth); Z[-1] = -100000

            for var, zgrid in zip(['temp', 'salt', 'u', 'v'], [zr, zr, zu, zv]):
                val = getattr(field, var)
                padded = np.vstack((val[0:1], val, val[-1:]))
                flipped = np.flip(padded, axis=0)
                sigma = tl.ztosigma_numba(flipped, zgrid, np.flipud(Z))
                setattr(field, var, sigma)

            # [12] Volume conservation
            print("--- [12] Volume conservation and barotropic correction ---")
            zw = tl.zlevs(*zargs, 5, grd.topo, field.zeta)
            dzr = zw[1:, :, :] - zw[:-1, :, :]
            dzu, dzv = tl.rho2uv(dzr, 'u'), tl.rho2uv(dzr, 'v')

            u_new, v_new, ubar_new, vbar_new = tl.conserve_and_recompute_barotropic(
                field.u, field.v, field.ubar, field.vbar, dzu, dzv
            )

            setattr(field, 'u', u_new)
            setattr(field, 'v', v_new)
            setattr(field, 'ubar', ubar_new)
            setattr(field, 'vbar', vbar_new)

            # [13] slicing to NSEW
            for varname in bry_data:
                var = getattr(field, varname)
                for direction in ['west', 'east', 'south', 'north']:
                    sliced = tl.extract_bry(var, direction)
                    bry_data[varname][direction].append(sliced)

            # 상대 시간도 저장
            time_converted = tl.compute_relative_time(tval, ogcm.time_unit, cfg.time_ref)
            bry_time.append(time_converted)
            print(num2date(time_converted,cfg.time_ref))
            end=time.time()
            print(f'--- Time elapsed: {end-start:.3f}s ---')
# --- 모든 시간 처리 후 최종 정리 ---
for varname in bry_data:
    for direction in bry_data[varname]:
        bry_data[varname][direction] = np.stack(bry_data[varname][direction], axis=0)

bry_time = np.array(bry_time)


# --- [10] Write all remapped variables to ini.nc ---
print(f"--- [13] Write all remapped variables to {cfg.bryname} ---")

with Dataset(cfg.bryname, 'a') as nc:
    # 시간 변수들
#    for tname in ['bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']:
#        nc[tname][:] = bry_time

    # 필드 저장
    for varname in bry_data:
        for direction in bry_data[varname]:
            var_fullname = f"{varname}_{direction}"
            nc[var_fullname][:] = bry_data[varname][direction]
print("--- [DONE] ---")


print(f'--- Time elapsed: {time.time()-start1:.3f}s ---')












