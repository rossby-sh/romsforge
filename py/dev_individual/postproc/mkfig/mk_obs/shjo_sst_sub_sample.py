# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import datetime
from netCDF4 import num2date, Dataset
from scipy.spatial import cKDTree
from matplotlib.path import Path
import os
# import writeObsfile  # 반드시 PYTHONPATH에 포함되어야 함
import time

# 경로 설정
pthG = 'D:/shjo/ROMS_inputs/'
pth = 'D:/shjo/GLORYS/OSTIA/'

# ROMS grid 열기
ncG = Dataset(pthG + 'roms_grd_fennel_5km_smooth_v2.nc')
lon_rho = ncG['lon_rho'][:]
lat_rho = ncG['lat_rho'][:]

# Observation 파일 출력 경로
outputFile = 'D:/shjo/ROMS_inputs/obs/obs_SST_OSTIA_10km.nc'
if os.path.exists(outputFile):
    os.remove(outputFile)

# NetCDF 열기
ds = xr.open_mfdataset(pth + "*.nc", decode_times=False).sel(
    lon=slice(110,160), lat=slice(5,55))

# 공통 설정
roms_ref = datetime.datetime(2000, 1, 1)
firstIteration = True
USENETCDF4 = True
Nstate = 12
obs_flag = 6  # SST

# 폴리곤 필터링 함수
def get_polygon(lon_rho, lat_rho):
    lon_bd = np.concatenate([lon_rho[:,0], lon_rho[-1,:], lon_rho[::-1,-1], lon_rho[0,::-1]])
    lat_bd = np.concatenate([lat_rho[:,0], lat_rho[-1,:], lat_rho[::-1,-1], lat_rho[0,::-1]])
    return np.vstack([lon_bd, lat_bd])

def filter_inside_roms_polygon(obs_lon, obs_lat, lon_rho, lat_rho):
    polygon = get_polygon(lon_rho, lat_rho)
    path = Path(polygon.T)
    obs_points = np.column_stack((obs_lon, obs_lat))
    return path.contains_points(obs_points)

def map_obs_to_roms_grid(obs_lon, obs_lat, lon_rho, lat_rho):
    tree = cKDTree(np.column_stack((lon_rho.flatten(), lat_rho.flatten())))
    obs_points = np.column_stack((obs_lon, obs_lat))
    _, idxs = tree.query(obs_points)
    j,i=np.unravel_index(idxs, lon_rho.shape)
    return i,j

def writeData(outputFile,obs_lat,obs_lon,obs_value,Nobs,survey_time,obs_time,obs_Xgrid,obs_Ygrid,
                               firstIteration,lastIteration,
                               obs_flag,obs_type,obs_error,obs_Zgrid,obs_depth,obs_variance,
                               survey,is3d,Nstate,USENETCDF4):
   if USENETCDF4 is True:
      myZLIB=True
      myformat='NETCDF4'
   else:
      myZLIB=False
      myformat='NETCDF3_CLASSIC'

   if firstIteration is True:

      f1 = Dataset(outputFile, mode='w', format=myformat)
      f1.description="This is a obs file for SST"
      f1.history = 'Created ' + time.ctime(time.time())
      f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
      f1.type='NetCDF4 using program createMapNS.py'

      f1.options='Program requires: getCortad.py and writeObsfile.py'

      """ Define dimensions """
      f1.createDimension('one', 1)
      f1.createDimension('state_variable', Nstate)
      f1.createDimension('datum', None)

      v_spherical = f1.createVariable('spherical', 'c', ('one',),zlib=myZLIB)
      v_spherical.long_name = 'grid type logical switch'
      v_spherical.option_T  = "spherical"
      v_spherical.option_F  = "Cartesian"
      v_spherical[:]        = "T"

      v_obs_type = f1.createVariable('obs_type', 'i', ('datum',),zlib=myZLIB)
      v_obs_type.long_name = 'model state variable associated with observation'
      v_obs_type.opt_1 ='free-surface'
      v_obs_type.opt_2 ='vertically integrated u-momentum component';
      v_obs_type.opt_3 ='vertically integrated v-momentum component';
      v_obs_type.opt_4 ='u-momentum component'
      v_obs_type.opt_5 ='v-momentum component'
      v_obs_type.opt_6 ='potential temperature'
      v_obs_type.opt_7 ='salinity'
      v_obs_type[:]    = obs_type

      v_time = f1.createVariable('obs_time', 'd', ('datum',),zlib=myZLIB)
      v_time.long_name = 'Time of observation'
      v_time.units     = 'days'
      v_time.field     = 'time, scalar, series'
      v_time.calendar  = 'standard'
      v_time[:]        = obs_time

      v_obs_lon = f1.createVariable('obs_lon', 'd', ('datum',),zlib=myZLIB)
      v_obs_lon.long_name = 'Longitude of observation'
      v_obs_lon.units     = 'degrees_east'
      v_obs_lon.min       = -180
      v_obs_lon.max       = 180
      v_obs_lon[:]        = obs_lon

      v_obs_lat = f1.createVariable('obs_lat', 'd', ('datum',),zlib=myZLIB)
      v_obs_lat.long_name = 'Latitude of observation'
      v_obs_lat.units     = 'degrees_north'
      v_obs_lat.min       = -90
      v_obs_lat.max       = 90
      v_obs_lat[:]        = obs_lat

      v_obs_depth = f1.createVariable('obs_depth', 'd', ('datum',),zlib=myZLIB)
      v_obs_depth.long_name = 'Depth of observation'
      v_obs_depth.units     = 'meter'
      v_obs_depth.minus     = 'downwards'
      v_obs_depth[:]        = obs_depth

      v_obs_error = f1.createVariable('obs_error', 'd', ('datum',),zlib=myZLIB)
      v_obs_error.long_name = 'Observation error covariance'
      v_obs_error.units     = 'squared state variable units'
      v_obs_error[:]        = obs_error

      v_obs_val = f1.createVariable('obs_value', 'd', ('datum',),zlib=myZLIB)
      v_obs_val.long_name = 'Observation value'
      v_obs_val.units     = 'state variable units'
      v_obs_val[:]        = obs_value

      v_obs_xgrid = f1.createVariable('obs_Xgrid', 'd', ('datum',),zlib=myZLIB)
      v_obs_xgrid.long_name = 'x-grid observation location'
      v_obs_xgrid.units     = 'nondimensional'
      v_obs_xgrid.left      = "INT(obs_Xgrid(datum))"
      v_obs_xgrid.right     = "INT(obs_Xgrid(datum))+1"
      v_obs_xgrid[:]        = obs_Xgrid

      v_obs_ygrid = f1.createVariable('obs_Ygrid', 'd', ('datum',),zlib=myZLIB)
      v_obs_ygrid.long_name = 'y-grid observation location'
      v_obs_ygrid.units     = 'nondimensional'
      v_obs_ygrid.top       = "INT(obs_Ygrid(datum))+1"
      v_obs_ygrid.bottom    = "INT(obs_Ygrid(datum))"
      v_obs_ygrid[:]        = obs_Ygrid

      v_obs_zgrid = f1.createVariable('obs_Zgrid', 'd', ('datum',),zlib=myZLIB)
      v_obs_zgrid.long_name = 'z-grid observation location'
      v_obs_zgrid.units     = 'nondimensional'
      v_obs_zgrid.up        = "INT(obs_Zgrid(datum))+1"
      v_obs_zgrid.down      = "INT(obs_Zgrid(datum))"
      v_obs_zgrid[:]        = obs_Zgrid

      t0 = time.time()
      """Find index for ading new info to arrays (same for all variables)"""
      myshape=f1.variables["obs_Zgrid"][:].shape
      indexStart=myshape[0]
      indexEnd=obs_Zgrid.shape[0]+myshape[0]
      t1 = time.time()
      print ("array append created in %s seconds"%(t1-t0))

      f1.close()

   # if firstIteration is False and lastIteration is False:
   if firstIteration is False :

      f1 = Dataset(outputFile, mode='a', format=myformat)

      t0 = time.time()
      """Find index for ading new info to arrays (same for all variables)"""
      myshape=f1.variables["obs_Zgrid"][:].shape
      indexStart=myshape[0]
      indexEnd=obs_Zgrid.shape[0]+myshape[0]

      f1.variables["obs_type"][indexStart:indexEnd] = obs_type
      f1.variables["obs_time"][indexStart:indexEnd] = obs_time
      f1.variables["obs_lon"][indexStart:indexEnd] = obs_lon
      f1.variables["obs_lat"][indexStart:indexEnd] = obs_lat
      f1.variables["obs_depth"][indexStart:indexEnd] = obs_depth
      f1.variables["obs_error"][indexStart:indexEnd] = obs_error
      f1.variables["obs_value"][indexStart:indexEnd] = obs_value
      f1.variables["obs_Xgrid"][indexStart:indexEnd] = obs_Xgrid
      f1.variables["obs_Ygrid"][indexStart:indexEnd] = obs_Ygrid
      f1.variables["obs_Zgrid"][indexStart:indexEnd] = obs_Zgrid
    
      t1 = time.time()
      print ("array append created in %s seconds"%(t1-t0))
      f1.close()

   if lastIteration is True:

      f1 = Dataset(outputFile, mode='a', format=myformat)

      f1.createDimension('survey', survey)

      v_obs = f1.createVariable('Nobs', 'i', ('survey',),zlib=myZLIB)
      v_obs.long_name = 'number of observations with the same survey time'
      v_obs.field     = 'scalar, series'
      v_obs[:]        = Nobs

      v_time = f1.createVariable('survey_time', 'i', ('survey',),zlib=myZLIB)
      v_time.long_name = 'Survey time'
      v_time.units     = 'day'
      v_time.field     = 'time, scalar, series'
      v_time.calendar  = 'standard'
      v_time[:]        = survey_time

      v_obs_var = f1.createVariable('obs_variance', 'd', ('state_variable',),zlib=myZLIB)
      v_obs_var.long_name = 'global time and space observation variance'
      v_obs_var.units     = 'squared state variable units'
      v_obs_var[:]        = obs_variance

      f1.close()

all_Nobs=[]
all_survey_times = []
# 시간 루프 시작
for t in range(len(ds.time)):

    print(t)
    sst = ds['analysed_sst'].isel(time=t)
    error = ds['analysis_error'].isel(time=t)
    lat = ds['lat']
    lon = ds['lon']

    # 유효값 마스크
    mask = (~np.isnan(sst)) & (~np.isnan(error)) & (sst < 400)
    sst = sst.where(mask)
    error = error.where(mask)

    lon2d, lat2d = np.meshgrid(lon, lat)

    obs_value = (sst - 273.15).values[::3,::3].flatten()
    obs_error = error.values[::3,::3].flatten()
    obs_lat = lat2d[::3,::3].flatten()
    obs_lon = lon2d[::3,::3].flatten()

    valid = (~np.isnan(obs_value)) & (~np.isnan(obs_error))
    obs_value = obs_value[valid]
    obs_error = obs_error[valid]
    obs_lat = obs_lat[valid]
    obs_lon = obs_lon[valid]

    # 시간 계산
    obs_time_raw = ds['time'].values[t]
    time_units = ds['time'].attrs['units']
    obs_datetime = num2date(obs_time_raw, units=time_units)
    obs_time = (obs_datetime - roms_ref).total_seconds() / 86400.0
    all_survey_times.append(obs_time)  # float
    obs_time = np.ones(len(obs_value)) * obs_time  # 관측마다 동일한 시간

    # 필터링
    inside = filter_inside_roms_polygon(obs_lon, obs_lat, lon_rho, lat_rho)
    obs_lon = obs_lon[inside]
    obs_lat = obs_lat[inside]
    obs_value = obs_value[inside]
    obs_error = obs_error[inside]
    obs_time = obs_time[inside]

    if len(obs_value) == 0:
        continue

    obs_Xgrid, obs_Ygrid = map_obs_to_roms_grid(obs_lon, obs_lat, lon_rho, lat_rho)

    Nobs = len(obs_value)
    all_Nobs.append(Nobs)
    obs_type = np.ones(Nobs, dtype=int) * obs_flag
    obs_depth = np.zeros(Nobs)
    obs_Zgrid = np.zeros(Nobs)

    survey_time = np.array([obs_time[0]])
    survey = t+1
    is3d = 0
    obs_variance = np.ones(Nstate)

    lastIteration = (t == len(ds.time) - 1)

    writeData(
        outputFile,
        obs_lat, obs_lon, obs_value, all_Nobs,
        all_survey_times, obs_time,
        obs_Xgrid, obs_Ygrid,
        firstIteration, lastIteration,
        obs_flag, obs_type, obs_error,
        obs_Zgrid, obs_depth, obs_variance,
        survey, is3d, Nstate, USENETCDF4
    )

    firstIteration = False
