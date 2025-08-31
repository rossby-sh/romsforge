# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 15:20:17 2025

@author: ust21
"""

# -*- coding: utf-8 -*-
import sys
import xarray as xr
import numpy as np
from datetime import timedelta
import datetime
from netCDF4 import num2date, Dataset
from scipy.spatial import cKDTree
from matplotlib.path import Path
import os
# import writeObsfile  # 반드시 PYTHONPATH에 포함되어야 함
import time
import pandas as pd
import matplotlib.pyplot as plt
PKG_path = 'C:/Users/ust21/shjo/projects/myROMS/prc_src/utils/' # Location of JNUROMS directory
sys.path.append(PKG_path)
from ROMS_utils01 import zlevs

# 경로 설정
pthG = 'D:/shjo/ROMS_inputs/NWP4_grd_3_10m_LP.nc'
pth = 'C:/Users/ust21/shjo/projects/myROMS/obs_data/'

# ROMS grid 열기
ncG = Dataset(pthG)
lon_rho = ncG['lon_rho'][:]
lat_rho = ncG['lat_rho'][:]

# Observation 파일 출력 경로
outputFile = 'D:/shjo/ROMS_inputs/TEST_obs.nc'
if os.path.exists(outputFile):
    os.remove(outputFile)
    
Table00=pd.read_excel(pth+'KODC_2023_2024.xlsx')
Table00.columns=['Region','_','_','_','Lat','Lon','Time','_','Depth','Temp','Sal','_','_','_','_','_','_','_','_','_']


Table00['Time'] = pd.to_datetime(Table00['Time'])


Table= Table00[['Lat','Lon','Depth','Temp','Sal']]

# 원하는 정보 추출
Table['Year'] = Table00['Time'].dt.year
Table['Month'] = Table00['Time'].dt.month
Table['Day'] = Table00['Time'].dt.day
Table['hour'] = Table00['Time'].dt.hour  # 'min'은 예약어라 이름을 다르게 쓰는 게 좋음
Table['min'] = Table00['Time'].dt.minute  # 'min'은 예약어라 이름을 다르게 쓰는 게 좋음



# 공통 설정
roms_ref = datetime.datetime(2000, 1, 1)
firstIteration = True
USENETCDF4 = True
Nstate = 7
all_Nobs=[]
all_survey_times = []
survey = 0
obs_flag = 999 
# 예시 파라미터 
Vtransform = 2
Vstretching = 2
theta_s = 7
theta_b = 0.1
hc = 200
N = 20
igrid = 1  # rho grid

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
        f1.description="This is a obs file for KODC"
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
        v_obs_prov.flag_values = np.array([1,2,3,4,5,6,7,8,9,10,11,12], dtype='float64')
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

        v_time = f1.createVariable('survey_time', 'f4', ('survey',),zlib=myZLIB)
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



# datetime.datetime(year, month, day, hour, minute)
datetimes = np.array([
    datetime.datetime(y, m, d, t) #- timedelta(hours=9) 
    for y, m, d, t in zip(Table['Year'].values, Table['Month'].values, Table['Day'].values, Table['hour'].values,)
])
Table['time']=datetimes

# 3. ROMS 기준 날짜로부터 days 계산
roms_ref = datetime.datetime(2000, 1, 1)
obs_time = np.array([(dt - roms_ref).total_seconds() / 86400.0 for dt in datetimes])


A=pd.DataFrame({'lat':Table['Lat'].values,'lon':Table['Lon'].values,'temp':Table['Temp'].values,'salt':Table['Sal'].values,'depth':Table['Depth'].values,'time':Table.time.values})
A['obs_time']=obs_time
df_sorted = A.sort_values(by="time")
df_sorted=df_sorted.set_index(df_sorted.time).drop('time',axis=1)
idt=np.unique(df_sorted.index)
#A=A.set_index(A.time).drop('time',axis=1)

for ii in df_sorted.groupby(level=0):
    df_t = ii[1]  # 시간 그룹의 DataFrame
    
    obs_lat = df_t['lat'].values
    obs_lon = df_t['lon'].values

    # 바운더리 필터
    inside = filter_inside_roms_polygon(obs_lon, obs_lat, lon_rho, lat_rho)
    df_inside = df_t[inside]

    # temp와 salt 모두 유효한 행들 추출 → 병합
    rows = []

    for _, row in df_inside.iterrows():
        if not np.isnan(row['temp']):
            rows.append({
                "lat": row['lat'],
                "lon": row['lon'],
                "depth": row['depth'],
                "value": row['temp'],
                "type": 6,
                "error": 0.3,
                "provenance": 4  # CTD temperature from Met Office
            })
        if not np.isnan(row['salt']):
            rows.append({
                "lat": row['lat'],
                "lon": row['lon'],
                "depth": row['depth'],
                "value": row['salt'],
                "type": 7,
                "error": 0.05,
                "provenance": 5  # CTD salinity from Met Office
            })

    if len(rows) == 0:
        continue

    # DataFrame으로 병합
    df_obs = pd.DataFrame(rows)

    # 값 추출
    obs_lat   = df_obs["lat"].values
    obs_lon   = df_obs["lon"].values
    obs_depth = df_obs["depth"].values
    obs_value = df_obs["value"].values
    obs_type  = df_obs["type"].values
    obs_provenance = df_obs["provenance"].values.astype("f8")
    # obs_error = np.maximum(df_obs["error"].values * np.abs(obs_value), 1e-3)  # 최소값 제한
    obs_error = np.maximum(df_obs["error"].values**2, 1e-5)  # 최소값 제한

    # obs_Zgrid = np.zeros_like(obs_value)

    # grid 매핑
    obs_Xgrid, obs_Ygrid = map_obs_to_roms_grid(obs_lon, obs_lat, lon_rho, lat_rho)

    # 시간 처리
    obs_time_val = round(df_inside['obs_time'].values[0],2)
    obs_time     = np.ones_like(obs_value) * obs_time_val


    ncG = Dataset(pthG)
    h=ncG['h'][:]
    h_point=h[obs_Ygrid, obs_Xgrid]
    
    # if not np.unique(h_point).shape:
    #     raise

    h_point_2d = h[obs_Ygrid, obs_Xgrid]
    zeta_2d = np.zeros_like(h_point_2d)
    z_r = zlevs(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h_point_2d, zeta_2d)
    # z_r: shape (N, M) → transpose to (M, N) if needed
    z_r = np.squeeze(z_r.T)  # shape (M, N)

    obs_Zgrid = np.array([
        int(np.argmin(np.abs(z_profile - z_obs)))
        for z_profile, z_obs in zip(z_r, -obs_depth)
    ])
    
    Nobs = len(obs_value)
    all_Nobs.append(Nobs)
    all_survey_times.append(obs_time_val)
    survey += 1
    survey_time = np.array([obs_time[0]])
    obs_variance =[0., 0., 0., 0., 0., 0.09, 0.0025 ]

    is3d = 1
    lastIteration = (survey == len(idt))
    
    writeData(
        outputFile,
        obs_lat, obs_lon, obs_value, all_Nobs,
        all_survey_times, obs_time,
        obs_Xgrid, obs_Ygrid,
        firstIteration, lastIteration,
        obs_flag, obs_type, obs_error,
        obs_Zgrid, obs_depth, obs_variance,obs_provenance,
        survey, is3d, Nstate, USENETCDF4
    )
    firstIteration = False
    
    

