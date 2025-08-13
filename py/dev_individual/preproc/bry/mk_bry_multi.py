

# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_B as cn
import utils as tl
from io_utils import collect_time_info
from collections import defaultdict
import time 

# --- [01] Load configuration and input metadata ---
cfg  = tl.parse_config("./config_multi.yaml")
grd  = tl.load_roms_grid(cfg.grdname)

filename_dict = {}
time_ref_list = None
tinfo_dict={}

for varname, meta in cfg.ogcm_inputs.to_dict().items():
    path = meta["path"]
    var = meta["varname"]
    filelist = sorted(glob.glob(os.path.join(path, "*.nc")))
    assert len(filelist) > 0, f"No files found for variable: {varname}"

    tinfo = collect_time_info(filelist, cfg.ogcm_var_name["time"], 
                               (str(cfg.bry_start_date), str(cfg.bry_end_date)))
    tinfo_dict[varname] = tinfo
    filename_dict[varname] = [t.filename for t in tinfo]

    # 공통 시간 비교
    current_dt_list = [t.datetime for t in tinfo]
    if time_ref_list is None:
        time_ref_list = current_dt_list
    else:
        if time_ref_list != current_dt_list:
            raise ValueError(f"[ERROR] Time mismatch in variable {varname}")


#for varname, flist in filename_dict.items():
#    print(f"[CHECK] {varname}: {len(flist)} files")
#    for f in flist[:3]:
#        print(f"  ↳ {f}")

for varname, tlist in tinfo_dict.items():
    print(f"[CHECK] {varname} tinfo count: {len(tlist)}")
    print("  First 3 time entries:")
    for t in tlist[:3]:
        print(f"    {t.datetime.isoformat()} from {os.path.basename(t.filename)}")


print("[CHECK] Time consistency across variables:")
for varname, tlist in tinfo_dict.items():
    mismatch = [t.datetime != ref for t, ref in zip(tlist, time_ref_list)]
    if any(mismatch):
        print(f"  [!ERROR] Time mismatch in {varname}")
        raise
    else:
        print(f"  [NOTE] {varname} aligned")

ogcm = tl.load_ogcm_metadata(filename_dict['zeta'][0], cfg.ogcm_var_name)
tinfo = tinfo_dict['zeta']

datenums = np.array([ti.raw_value for ti in tinfo])
relative_time = tl.compute_relative_time(datenums, ogcm.time_unit, cfg.time_ref)

# --- [03] Create initial NetCDF file ---
print(f"--- [03] Creating boundary NetCDF file: {cfg.bryname} ---")
status = cn.create_bry(cfg, grd, relative_time,  bio_model=None, ncFormat=cfg.ncformat)
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


bry_data = tl.make_all_bry_data_shapes(['zeta','ubar','vbar','temp','salt','u','v'], len(tinfo), grd, cfg.vertical.layer_n)

bry_time = []

# --- [05] Load OGCM raw fields (zeta, temp, salt, u, v) ---
print("--- [05] Loading OGCM data ---")
var_list = ['zeta', 'temp', 'salt', 'u', 'v']
filename_tuple_list = list(zip(*[filename_dict[v] for v in var_list]))
print(filename_tuple_list)

grouped = defaultdict(list)
for entry in tinfo:
    grouped[entry.filename].append(entry)

# --- 정렬 추가 ---
for entries in grouped.values():
    entries.sort(key=lambda x: x.datetime)

time_index_map = {entry.datetime: n for n, entry in enumerate(tinfo)}
import time

# --- 파일 단위 루프 ---
for (f_zeta, f_temp, f_salt, f_u, f_v), entries in zip(filename_tuple_list, grouped.values()):
    nc_zeta = tl.MaskedNetCDF(Dataset(f_zeta, maskandscale=True))
    nc_temp = tl.MaskedNetCDF(Dataset(f_temp, maskandscale=True))
    nc_salt = tl.MaskedNetCDF(Dataset(f_salt, maskandscale=True))
    nc_u    = tl.MaskedNetCDF(Dataset(f_u,    maskandscale=True))
    nc_v    = tl.MaskedNetCDF(Dataset(f_v,    maskandscale=True))

    for entry in entries:
        i = entry.index
        t = entry.datetime
        tval = entry.raw_value
        n = time_index_map[t]

        zeta = nc_zeta.get(cfg.ogcm_inputs["zeta"]['varname'], i, idy, idx)
        temp = nc_temp.get(cfg.ogcm_inputs["temp"]['varname'], i, slice(None), idy, idx)
        salt = nc_salt.get(cfg.ogcm_inputs["salt"]['varname'], i, slice(None), idy, idx)
        u    = nc_u.get(cfg.ogcm_inputs["u"]['varname'], i, slice(None), idy, idx)
        v    = nc_v.get(cfg.ogcm_inputs["v"]['varname'], i, slice(None), idy, idx)
        ubar = tl.depth_average(u, ogcm.depth)
        vbar = tl.depth_average(v, ogcm.depth)

        field = tl.ConfigObject(zeta=zeta, ubar=ubar, vbar=vbar,
                                temp=temp, salt=salt, u=u, v=v)
        start=time.time()
        # [06] Remap
        print("--- [06] Remapping ---")
        for var in vars(field):
            var_src = getattr(field, var)
            remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
            setattr(field, var, remapped)

        print(f'--- Time elapsed: {time.time()-start:.3f}s ---')
        start=time.time()
        
        # [07] Horizontal Flood
        print("--- [07] Apply horizontal flood ---")
        for var in vars(field):
            val = getattr(field, var)
            val_flooded = tl.flood_horizontal(val, grd.lon, grd.lat, method=cfg.flood_method_for_bry)
            setattr(field, var, val_flooded)
        print(f'--- Time elapsed: {time.time()-start:.3f}s ---')
        start=time.time()

        print("--- [08] Apply vertical flood ---")

        for var in vars(field):
            val = getattr(field, var)
            if val.ndim==2:
                continue
            #val_flooded = tl.flood_vertical_vectorized(val, grd.mask, spval=-1e10)
            val_flooded = tl.flood_vertical_numba(np.asarray(val), np.asarray(grd.mask), spval=-1e10)
            setattr(field, var, val_flooded)
        print(f'--- Time elapsed: {time.time()-start:.3f}s ---')
        start=time.time()
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
        print(f'--- Time elapsed: {time.time()-start:.3f}s ---')
        start=time.time()

        # [11] Vertical interpolation
        print("--- [11] Vertical interpolation from z to sigma ---")

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
        dz_3d = zw_3d[1:, :, :] - zw_3d[:-1, :, :]

        zu_3d = tl.rho2uv(zr_3d,"u")
        zv_3d = tl.rho2uv(zr_3d,"v")
        dz_u3d = tl.rho2uv(dz_3d,"u")
        dz_v3d = tl.rho2uv(dz_3d,"v")

        for direction in directions:
            zr = tl.extract_bry(zr_3d, direction)
            zu = tl.extract_bry(zu_3d, direction)
            zv = tl.extract_bry(zv_3d, direction)
            dzr = tl.extract_bry(dz_3d, direction)
            dzu = tl.extract_bry(dz_u3d, direction)
            dzv = tl.extract_bry(dz_v3d, direction)

# 이후 loop에서 var에 따라 zr/zu/zv 선택해서 사용


            for var, zgrid in zip(['temp', 'salt', 'u', 'v'],[zr,zr,zu,zv]):
                val = tl.extract_bry(getattr(field, var), direction)        # (Nz, Lp)
                val_pad = np.vstack((val[0:1], val, val[-1:]))              # (Nz+2, Lp)
                val_flip = np.flip(val_pad, axis=0)                         # (Nz+2 → 위에서 아래로)
                sigma_interped = tl.ztosigma_1d_numba(val_flip, zgrid, Z_flipped)             # (Ns, Lp)
                
                if var == 'u':
                    barotropic = tl.extract_bry(field.ubar, direction)
                    dz_bar = dzu
                elif var == 'v':
                    barotropic = tl.extract_bry(field.vbar, direction)
                    dz_bar = dzv
                else:
                    barotropic = None

                if barotropic is not None:
                    sigma_interped, _, barotropic_corrected, _ = tl.conserve_and_recompute_barotropic(
                            sigma_interped, sigma_interped, barotropic, barotropic, dz_bar, dz_bar)
                    setattr(field, f"{'ubar' if var == 'u' else 'vbar'}_{direction}", barotropic_corrected)

                setattr(field, f"{var}_{direction}", sigma_interped)

            # zeta만 따로 저장 (ubar/vbar는 위에서 이미 저장됨)
            val = tl.extract_bry(field.zeta, direction)
            setattr(field, f"zeta_{direction}", val)

            #Save to data
            for varname in bry_data:
                val_d = getattr(field, f"{varname}_{direction}")
                bry_data[varname][direction][n, ...] = val_d

            # 시간 저장
            time_converted = tl.compute_relative_time(tval, ogcm.time_unit, cfg.time_ref)
            bry_time.append(time_converted)
        print(num2date(time_converted,cfg.time_ref)) 

    nc_zeta.nc.close()
    nc_temp.nc.close()
    nc_salt.nc.close()
    nc_u.nc.close()
    nc_v.nc.close()

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














