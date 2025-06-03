import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
#from py.dev_individual.libs.tools import crop_to_model_domain, load_ogcm_metadata
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import create_I as cn
import tools as tl
import xesmf as xe
import xarray as xr
from io_utils import collect_time_info

cfg = tl.parse_config("./config.yaml")

grd = tl.load_roms_grid(cfg.grdname)
ogcm = tl.load_ogcm_metadata(cfg.ogcm_name, cfg.ogcm_var_name)


tinfo = collect_time_info(cfg.ogcm_name, cfg.ogcm_var_name['time'], cfg.initdate, units=cfg.time_ref)
_, time_index, _ = tinfo[0]

relative_time = tl.compute_relative_time(ogcm.time[time_index], ogcm.time_unit, cfg.time_ref)

status = cn.create_ini(cfg, grd, relative_time, ncFormat=cfg.ncformat, bio_model="Fennel")
if status == 0:
    print(f'--- Failed to creating file {cfg.ininame}')
    raise


lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)

if cfg.calc_weight:
    print(f'--- Calculating weight: {cfg.weight_file} ---')
    status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
    
    if status == 0:
        print("☠️리매핑 weight 생성 실패 → 스킵 또는 에러 로그")
        raise  # or raise
else: 
    print(f'--- Use existing wght file {cfg.weight_file} ---')


'''
## Read data and and weight components
 RUTERGUS tools 처럼 변수 차원 확인해서 자동으로 슬라이싱 
'''


# 필요한 정보 추출
file_path, time_index, _ = tinfo[0]

# 데이터 슬라이싱
with Dataset(file_path) as nc:
    zeta = nc[cfg.ogcm_var_name['zeta']][time_index, idy, idx]
    temp = nc[cfg.ogcm_var_name['temperature']][time_index, :, idy, idx]
    salt = nc[cfg.ogcm_var_name['salinity']][time_index, :, idy, idx]
    u    = nc[cfg.ogcm_var_name['u']][time_index, :, idy, idx]
    v    = nc[cfg.ogcm_var_name['v']][time_index, :, idy, idx]


ubar, vbar 구하는 함수 추가 
변수들 객체화하여 다루기







with Dataset(cfg.weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]




'''
Apply weight based interpolation
'''
print(zeta.shape)
var_src=zeta.data
ydim, xdim = grd.lon.shape
assert var_src.ndim == 2
src_flat = var_src.flatten()
dst_flat = np.zeros(row.max() + 1)

np.add.at(dst_flat, row, S * src_flat[col])

dst_flat=dst_flat.reshape( (ydim, xdim) )

# debugging 

print('--- Remapping temp ---')

var_src = temp.data
# shape 체크 (안전장치)
assert var_src.ndim == 3
nz, ny, nx = var_src.shape
src_flat = var_src.reshape(nz, ny * nx)
var_dst=np.zeros([nz,ydim, xdim])

# remapping 루프
for k in range(nz):
    dst_flat = np.zeros(row.max() + 1)
    np.add.at(dst_flat, row, S * src_flat[k, col])
    var_dst[k] = dst_flat.reshape( (ydim, xdim) )












