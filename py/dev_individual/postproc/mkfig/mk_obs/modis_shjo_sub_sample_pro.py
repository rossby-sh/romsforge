# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 09:00:00 2025
@author: ust21
"""
import glob
import xarray as xr
import numpy as np
import datetime
from netCDF4 import Dataset, num2date
from scipy.spatial import cKDTree
from matplotlib.path import Path
import os
import time
# import writeObsfile  # PYTHONPATH에 포함 필요

# 경로 설정
pthG = 'D:/shjo/ROMS_inputs/'
pth = 'D:/shjo/MODIS/AQUA_chl_NRT/'
outputFile = 'D:/shjo/ROMS_inputs/obs/pro/obs_phyt_27km.nc'

# ROMS grid 열기
ncG = Dataset(pthG + 'roms_grd_fennel_15km_smooth_v2.nc')
lon_rho = ncG['lon_rho'][:]
lat_rho = ncG['lat_rho'][:]


files = sorted(glob.glob(pth + '*.nc'))
datasets = []
ref_date = np.datetime64('2000-01-01')

for f in files:
    ds = xr.open_dataset(f)
    ds = ds.expand_dims('time')

    # 날짜 추출
    fname = os.path.basename(f)
    date_str = fname.split('.')[1]  # 예: '20241231'
    dt = datetime.datetime.strptime(date_str, "%Y%m%d")
    np_time = (np.datetime64(dt) - ref_date) / np.timedelta64(1, 'D')  # float days

    ds['time'] = ('time', [np_time])
    ds['time'].attrs['units'] = 'days since 2000-01-01 00:00:00'
    datasets.append(ds)

ds_combined = xr.concat(datasets, dim='time').sel(
    lon=slice(100,170), lat=slice(55,5))

# 공통 설정
roms_ref = datetime.datetime(2000, 1, 1)
firstIteration = True
USENETCDF4 = True
Nstate = 19
obs_flag = 11  # CHL-a
obs_provenance=14

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
               obs_flag,obs_type,obs_error,obs_Zgrid,obs_depth,obs_variance,obs_provenance,
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
        v_obs_type.opt_1  = 'free-surface'
        v_obs_type.opt_2  = 'vertically integrated u-momentum component'
        v_obs_type.opt_3  = 'vertically integrated v-momentum component'
        v_obs_type.opt_4  = 'u-momentum component'
        v_obs_type.opt_5  = 'v-momentum component'
        v_obs_type.opt_6  = 'potential temperature'
        v_obs_type.opt_7  = 'salinity'
        v_obs_type.opt_8  = 'NH4'
        v_obs_type.opt_9  = 'NO3'
        v_obs_type.opt_10 = 'chlorophyll'
        v_obs_type.opt_11 = 'phytoplankton'
        v_obs_type.opt_12 = 'zooplankton'
        v_obs_type.opt_13 = 'LdetritusN'
        v_obs_type.opt_14 = 'SdetritusN'
        v_obs_type.opt_15 = 'oxygen'
        v_obs_type.opt_16 = 'PO4'
        v_obs_type.opt_17 = 'LdetritusP'
        v_obs_type.opt_18 = 'SdetritusP'
        v_obs_type.opt_19 = 'H2S'
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

        v_obs_prov = f1.createVariable('obs_provenance', 'd', ('datum',), zlib=myZLIB)
        v_obs_prov.long_name = 'observation origin'
        v_obs_prov.flag_values = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14], dtype='float64')
        v_obs_prov.opt_1  = "gridded AVISO sea level anomaly"
        v_obs_prov.opt_2  = "blended satellite SST"
        v_obs_prov.opt_3  = "XBT temperature from Met Office"
        v_obs_prov.opt_4  = "CTD temperature from Met Office"
        v_obs_prov.opt_5  = "CTD salinity from Met Office"
        v_obs_prov.opt_6  = "ARGO floats temperature from Met Office"
        v_obs_prov.opt_7  = "ARGO floats salinity from Met Office"
        v_obs_prov.opt_8  = "CTD temperature from CalCOFI"
        v_obs_prov.opt_9  = "CTD salinity from CalCOFI"
        v_obs_prov.opt_10 = "CTD temperature from GLOBEC"
        v_obs_prov.opt_11 = "CTD salinity from GLOBEC"
        v_obs_prov.opt_12 = "buoy, thermistor temperature from Met Office"
        v_obs_prov.opt_13 = "CTD Chlorophyll from NIFS"
        v_obs_prov.opt_14 = "Satellite Chlorophyll"
        v_obs_prov[:] = obs_provenance

        f1.close()

    if firstIteration is False:
        f1 = Dataset(outputFile, mode='a', format=myformat)

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
        f1.variables["obs_provenance"][indexStart:indexEnd] = obs_provenance

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




all_survey_times = []
all_Nobs=[]

if os.path.exists(outputFile):
    os.remove(outputFile)

for t in range(len(ds_combined.time)):
#for t in range(5):

    print("Processing time index:", t)
    chl = ds_combined['chlor_a'].isel(time=t)
    lat = ds_combined['lat']
    lon = ds_combined['lon']

    # 유효값 필터링
    mask = (chl > 0.01) & (chl < 100) & ~np.isnan(chl)
    chl = chl.where(mask)

    lon2d, lat2d = np.meshgrid(lon, lat)

    obs_value = chl.values[::3,::3].flatten()/(0.02*6.625*12)  # 단위 변환: mg/m^3 → mmol/m^3
    obs_lat = lat2d[::3,::3].flatten()
    obs_lon = lon2d[::3,::3].flatten()

    valid = ~np.isnan(obs_value)
    obs_value = obs_value[valid]
    obs_lat = obs_lat[valid]
    obs_lon = obs_lon[valid]

    # 시간 계산
    obs_time_raw = ds_combined['time'].values[t]
    time_units = ds_combined['time'].attrs['units']
    obs_datetime = num2date(obs_time_raw, units=time_units)
    obs_time_val = (obs_datetime - roms_ref).total_seconds() / 86400.0
    all_survey_times.append(obs_time_val)
    obs_time = np.ones(len(obs_value)) * obs_time_val
    
    # 도메인 필터링
    inside = filter_inside_roms_polygon(obs_lon, obs_lat, lon_rho, lat_rho)
    obs_lon = obs_lon[inside]
    obs_lat = obs_lat[inside]
    obs_value = obs_value[inside]
    obs_time = obs_time[inside]

    if len(obs_value) == 0:
        continue
    

    obs_Xgrid, obs_Ygrid = map_obs_to_roms_grid(obs_lon, obs_lat, lon_rho, lat_rho)
    
    Nobs = len(obs_value)
    all_Nobs.append(Nobs)
    obs_type = np.ones(Nobs, dtype=int) * obs_flag
    obs_error = np.maximum(obs_value * 0.03, 1e-5)
    obs_depth = np.zeros(Nobs)
    obs_Zgrid = np.zeros(Nobs)
    is3d = 0
    obs_variance = np.ones(Nstate)

    survey = t+1
    lastIteration = (t == len(ds_combined.time) - 1)

    survey_time = np.array(all_survey_times) if lastIteration else np.array([obs_time[0]])
    
    writeData(
        outputFile,
        obs_lat, obs_lon, obs_value, all_Nobs,
        survey_time, obs_time,
        obs_Xgrid, obs_Ygrid,
        firstIteration, lastIteration,
        obs_flag, obs_type, obs_error,
        obs_Zgrid, obs_depth, obs_variance,obs_provenance,
        survey, is3d, Nstate, USENETCDF4
    )

    firstIteration = False
    
    
    
