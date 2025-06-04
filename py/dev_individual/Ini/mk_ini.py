import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
#from py.dev_individual.libs.tools import crop_to_model_domain, load_ogcm_metadata
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import create_I as cn
import utils as tl
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
if status:
    print(f'--- Failed to creating file {cfg.ininame}')
    raise


lon_crop, lat_crop, idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)

if cfg.calc_weight:
    print(f'--- Calculating weight: {cfg.weight_file} ---')
    status = tl.build_bilinear_regridder(lon_crop, lat_crop, grd.lon, grd.lat, cfg.weight_file, reuse=False)
    
    if status:
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



with Dataset(file_path, maskandscale=False) as nc_raw:
    nc = tl.MaskedNetCDF(nc_raw)

    zeta = nc.get(cfg.ogcm_var_name['zeta'], time_index, idy, idx)
    temp = nc.get(cfg.ogcm_var_name['temperature'], time_index, slice(None), idy, idx)
    salt = nc.get(cfg.ogcm_var_name['salinity'], time_index, slice(None), idy, idx)
    u    = nc.get(cfg.ogcm_var_name['u'], time_index, slice(None), idy, idx)
    v    = nc.get(cfg.ogcm_var_name['v'], time_index, slice(None), idy, idx)

    ubar = tl.depth_average(u, ogcm.depth)
    vbar = tl.depth_average(v, ogcm.depth)

ogcm_data = tl.ConfigObject(
    zeta = zeta,
    ubar = ubar,
    vbar = vbar,
    temp = temp,
    salt = salt,
    u    = u,
    v    = v
)

with Dataset(cfg.weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]


for varname in vars(ogcm_data) :
    var_src = getattr(ogcm_data, varname)
    remapped = tl.remap_variable(var_src, row, col, S, grd.lon.shape, method="coo")
    setattr(ogcm_data, varname, remapped)


# --- Vector processing ---

u=ru2.rho2u_3d(OGCM_Data['u']*cosa+OGCM_Data['v']*sina)
v=ru2.rho2v_3d(OGCM_Data['v']*cosa-OGCM_Data['u']*sina)

# Euler method
ubar=ru2.rho2u_2d(OGCM_Data['ubar']*cosa+OGCM_Data['vbar']*sina)
vbar=ru2.rho2v_2d(OGCM_Data['vbar']*cosa-OGCM_Data['ubar']*sina)

# --- Process ROMS vertical grid ---




