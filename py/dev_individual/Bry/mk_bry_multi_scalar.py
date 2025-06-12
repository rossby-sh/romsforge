

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

# --- [01] Load configuration and input metadata ---
cfg  = tl.parse_config("./config_multi.yaml")
grd  = tl.load_roms_grid(cfg.grdname)

cfg.ogcm_scalar_path = {
    k: v for k, v in cfg.ogcm_path.items() if k not in ['u', 'v']
}

"""
status = cn.create_bry(cfg, grd, relative_time,  bio_model=cfg.bio_model_type, ncFormat=cfg.ncformat)
if status:
    print(f"--- [!ERROR] Failed to creating file {cfg.bryname} ---")
    raise
print(f"--- Created file: {cfg.bryname} ---")
"""

for varname, path in cfg.ogcm_path.items():
    print(f"--- [VAR] {varname} from {path}")

    filelist = sorted(glob.glob(os.path.join(cfg.ogcm_path, "*.nc")))

    ogcm = tl.load_ogcm_metadata(filelist[0], cfg.ogcm_var_name)

    tinfo = collect_time_info(filelist, cfg.ogcm_var_name.time,\
            (str(cfg.bry_start_date), str(cfg.bry_end_date)) )

    datenums = np.array([tval for _, _, _, tval in tinfo])

    relative_time = tl.compute_relative_time(datenums, ogcm.time_unit, cfg.time_ref)

    lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)

    if cfg.calc_weight:
        print(f"--- Calculating weight: {cfg.weight_file} ---")
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


    bry_data = tl.make_bry_data_shape(varname, num_steps, grd, Ns)
    """
        bry_data[varname] = {
        'west':  np.zeros((num_steps, Mp), dtype=np.float32),
        'east':  np.zeros((num_steps, Mp), dtype=np.float32),
        'south': np.zeros((num_steps, Lp), dtype=np.float32),
        'north': np.zeros((num_steps, Lp), dtype=np.float32),
    }
    이런 배열을 만들어 줌. single ,multi file 둘다 대응하게 만들기
    """
    
    bry_time = []

    grouped = defaultdict(list)
    for f, i, t, tval in tinfo:
        grouped[f].append((i, t, tval))

    for f in grouped:
        grouped[f].sort(key=lambda x: x[1])  # datetime 기준 정렬


    # --- 파일 단위 루프 ---
    for f, entries in grouped.items():
        
        with Dataset(f, maskandscale=True) as nc:
            nc_wrap = tl.MaskedNetCDF(nc)
            for i, t, tval in entries:
                
                var_nc = nc.variables[cfg.ogcm_var_name[varname]]
                ndim = var_nc.ndim

                if ndim == 3:
                    data = nc_wrap.get(cfg.ogcm_var_name[varname], i, slice(None), idy, idx)  # 3D
                elif ndim == 2:
                    data = nc_wrap.get(cfg.ogcm_var_name[varname], i, idy, idx)               # 2D
                else:
                    raise ValueError(f"Unsupported ndim={ndim} for variable '{varname}'")
                
                field = tl.ConfigObject(**{varname: data})

                print("--- Remapping ---")
                for varname in vars(field):
                    var_src = getattr(field, varname)
                    remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
                    setattr(field, varname, remapped)

                print("--- Apply horizontal flood ---")
                for var in vars(field):
                    val = getattr(field, var)
                    val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_bry)
                    setattr(field, var, val_flooded)

                print("--- Apply vertical flood ---")
                for var in vars(field):
                    val = getattr(field, var)
                    if len(val.shape)==2:
                        continue
                    #val_flooded = tl.flood_vertical_vectorized(val, grd.mask, spval=-1e10)
                    val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask), spval=-1e10)
                    setattr(field, var, val_flooded)
               
                print("--- Masking land to 0 ---")
                for var in vars(field):
                    arr = getattr(field, var)
                    if arr.ndim == 2:
                        arr[grd.mask == 0] = 0.0
                    elif arr.ndim == 3:
                        arr[:, grd.mask == 0] = 0.0
                    setattr(field, var, arr)

                # [11] Vertical interpolation
                print("--- Vertical interpolation from z to sigma ---")

                # --- [11~13] 방향별 수직 보간 및 저장 ---
                directions = ['north', 'south', 'east', 'west']

                # HYCOM depth padding
                Z = np.zeros(len(ogcm.depth) + 2)
                Z[0] = 100
                Z[1:-1] = -ogcm.depth
                Z[-1] = -100000
                Z_flipped = np.flipud(Z)

                # sigma 관련 파라미터
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

                for direction in directions:
                    zr = tl.extract_bry(zr_3d, direction)

                    for var in vars(field):
                        val = tl.extract_bry(getattr(field, var), direction)
                        val_pad = np.vstack((val[0:1], val, val[-1:]))
                        val_flip = np.flip(val_pad, axis=0)
                        sigma_interped = tl.ztosigma_1d_numba(val_flip, zr, Z_flipped)
                        setattr(field, f"{var}_{direction}", sigma_interped)






                    # zeta만 따로 저장 (ubar/vbar는 위에서 이미 저장됨)
                    val = tl.extract_bry(field.zeta, direction)
                    setattr(field, f"zeta_{direction}", val)

                    # 저장
                    for varname in bry_data:
                        val_d = getattr(field, f"{varname}_{direction}")
                        bry_data[varname][direction].append(val_d)

                    # 시간 저장
                    time_converted = tl.compute_relative_time(tval, ogcm.time_unit, cfg.time_ref)
                    bry_time.append(time_converted)



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












