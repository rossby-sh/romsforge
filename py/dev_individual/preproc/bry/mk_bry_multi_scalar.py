

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
time_var = cfg.time_reference_var
filelist = sorted(glob.glob(os.path.join(cfg.ogcm_inputs[time_var]["path"], "*.nc")))

ogcm = tl.load_ogcm_metadata(filelist[0], cfg.ogcm_var_name)

tinfo = collect_time_info(filelist, cfg.ogcm_var_name.time,\
        (str(cfg.bry_start_date), str(cfg.bry_end_date)))
datenums = np.array([tval for _, _, _, tval in tinfo])
relative_time = tl.compute_relative_time(datenums, ogcm.time_unit, cfg.time_ref)

status = cn.create_bry(cfg, grd, relative_time,  bio_model=None, ncFormat=cfg.ncformat)
if status:
    print(f"--- [!ERROR] Failed to creating file {cfg.bryname} ---")
    raise
print(f"--- Created file: {cfg.bryname} ---")

for varname, meta in cfg.ogcm_inputs.to_dict().items():
    if varname in ['u','v']:
        print('--- [NOTE] Skip vectors ---')
        continue

    path = meta["path"]
    ogcm_varname = meta["varname"]
    print(f"--- [VAR] {varname} from {path}")

    filelist = sorted(glob.glob(os.path.join(path, "*.nc")))

    ogcm = tl.load_ogcm_metadata(filelist[0], cfg.ogcm_var_name)

    tinfo = collect_time_info(filelist, cfg.ogcm_var_name.time,\
            (str(cfg.bry_start_date), str(cfg.bry_end_date)) )

    time_index_map = {t: n for n, (_, _, t, _) in enumerate(tinfo)}

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


    
    bry_data = tl.make_bry_data_shape(varname, len(tinfo), grd, cfg.vertical.layer_n)
   

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
                n = time_index_map[t]
                var_nc = nc.variables[ogcm_varname]
                ndim = var_nc.ndim
                if ((ndim == 3) and (var_nc.shape[0] != 1)) or ((ndim == 4) and (var_nc.shape[0] == 1)) :
                    data = nc_wrap.get(ogcm_varname, i, slice(None), idy, idx)  # 3D
                elif (ndim == 2) or ((ndim == 3) and (var_nc.shape[0] == 1 )) :
                    data = nc_wrap.get(ogcm_varname, i, idy, idx)               # 2D
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
                Z[1:-1] = -np.abs(ogcm.depth)
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

                zr_3d = tl.zlevs(*zargs, 1, grd.topo, np.zeros_like(grd.topo))    # (Ns, Mp, Lp)
                zw_3d = tl.zlevs(*zargs, 5, grd.topo, np.zeros_like(grd.topo))    # (Ns+1, Mp, Lp)

                for var in list(vars(field)):
                    arr = getattr(field, var)

                    for direction in directions:
                    
                        if arr.ndim == 2:
                            val = tl.extract_bry(arr, direction)
                            #setattr(field, f"{direction}", val)
                        else :

                            zr = tl.extract_bry(zr_3d, direction)
                            val = tl.extract_bry(arr, direction)
                            val_pad = np.vstack((val[0:1], val, val[-1:]))
                            val_flip = np.flip(val_pad, axis=0)
                            val = tl.ztosigma_1d_numba(val_flip, zr, Z_flipped)
                            #setattr(field, f"{direction}", sigma_interped)
            
                        #if hasattr(field, direction):
                        print(f"[WRITE] var={varname}, direction={direction}, time index n={n}")
                        bry_data[direction][n, ...] = val

            print(f"---------------- test {varname} ------------------------")


    with Dataset(cfg.bryname, 'a') as nc:
        for direction, values in bry_data.items():
            nc_var = f"{varname}_{direction}"  # varname은 루프 바깥에서 정의되어 있음
            print(f"[WRITE] {nc_var} shape={values.shape}")
            nc.variables[nc_var][:] = values






