# -*- coding: utf-8 -*-
"""
Created on Tue May 13 13:53:14 2025

@author: ust21
"""

from netCDF4 import Dataset
from scipy  import io
import numpy as np

pth='D:/shjo/ROMS_inputs/NWP12_NWP4/'
outputFile=pth+'obs_sst_50km_N36_XY.nc'

NC=Dataset(pth+'obs_sst_50km_N36.nc','r')

MAT=io.loadmat(pth+'XYgrid_sst_50km.mat')

Xgrid=MAT['Xgrid'][:].reshape(-1)
Ygrid=MAT['Ygrid'][:].reshape(-1)

NC.variables.keys()

obs_type=NC['obs_type'][:]
obs_time=NC['obs_time'][:]
obs_lon=NC['obs_lon'][:]
obs_lat=NC['obs_lat'][:]
obs_depth=NC['obs_depth'][:]
obs_error=NC['obs_error'][:]
obs_value=NC['obs_value'][:]
obs_Zgrid=NC['obs_Zgrid'][:]
obs_provenance=NC['obs_provenance'][:]

idna=Xgrid==Xgrid 

obs_type=obs_type[idna]
obs_time=obs_time[idna]
obs_lon=obs_lon[idna]
obs_lat=obs_lat[idna]
obs_depth=obs_depth[idna]
obs_error=obs_error[idna]
obs_value=obs_value[idna]
obs_Zgrid=obs_Zgrid[idna]
obs_Xgrid=Xgrid[idna]
obs_Ygrid=Ygrid[idna]
obs_provenance=obs_provenance[idna]
    
survey_time, counts = np.unique(obs_time, return_counts=True)



Nobs=counts
survey_time= np.unique(obs_time)
obs_variance=NC['obs_variance'][:]

# NC['obs_type'][:]=obs_type
# NC['obs_time'][:]=obs_time
# NC['obs_lon'][:]=obs_lon
# NC['obs_lat'][:]=obs_lat
# NC['obs_depth'][:]=obs_depth
# NC['obs_error'][:]=obs_error
# NC['obs_value'][:]=obs_value
# NC['obs_Zgrid'][:]=obs_Zgrid
# NC['obs_Xgrid'][:]=obs_Xgrid
# NC['obs_Ygrid'][:]=obs_Ygrid

print(NC.dimensions['survey'].size)  # 혹은 NC.dimensions['Nsurvey'].size
print(len(survey_time), len(counts)) # 이것과 일치해야 함

NC.close()

Nstate=19
myZLIB=True

import time


f1 = Dataset(outputFile, mode='w', format='NETCDF3_CLASSIC')
f1.description="This is a obs file for OSTIA SST"
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
v_obs_type.opt_10  = 'chlorophyll'
v_obs_type.opt_11  = 'phytoplankton'
v_obs_type.opt_12  = 'zooplankton'
v_obs_type.opt_13  = 'LdetritusN'
v_obs_type.opt_14  = 'SdetritusN'
v_obs_type.opt_15  = 'oxygen'
v_obs_type.opt_16  = 'PO4'
v_obs_type.opt_17  = 'LdetritusP'
v_obs_type.opt_18  = 'SdetritusP'
v_obs_type.opt_19  = 'H2S'


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


f1.createDimension('survey', len(Nobs))

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







def add_global_attrs_to_obs(ncfiles):
    if isinstance(ncfiles, str):
        ncfiles = [ncfiles]

    for ncfile in ncfiles:
        print(f"🔧 Updating: {ncfile}")
        with Dataset(ncfile, 'a') as ds:
            # ===== 글로벌 속성 =====
            ds.title = 'North Western Pacific 25 km'
            ds.grd_file = 'roms_grd_fennel_15km_smooth_v2'
            ds.grid_Lm_Mm_N = np.array([563, 494, 36], dtype='i4')
            ds.variance_units = 'squared state variable units'

            ds.state_variables = "\n".join([
                "1: free-surface (m)",
                "2: vertically integrated u-momentum component (m/s)",
                "3: vertically integrated v-momentum component (m/s)",
                "4: u-momentum component (m/s)",
                "5: v-momentum component (m/s)",
                "6: temperature (Celsius)",
                "7: salinity",
                "8: NH4",
                "9: NO3",
                "10: chlorophyll",
                "11: phytoplankton",
                "12: zooplankton",
                "13: LdetritusN",
                "14: SdetritusN",
                "15: oxygen",
                "16: PO4",
                "17: LdetritusP",
                "18: SdetritusP",
                "19: H2S"
                
            ])

            ds.obs_provenance = "\n".join([
                "1: AVISO SLA",
                "2: blended SST",
                "3: XBT temperature",
                "4: CTD temperature",
                "5: CTD salinity",
                "6: ARGO temperature",
                "7: ARGO salinity",
                "8: CalCOFI temperature",
                "9: CalCOFI salinity",
                "10: GLOBEC temperature",
                "11: GLOBEC salinity",
                "12: buoy temp",
                "13: CTD Chlorophyll",
                "14: Satellite Chlorophyll"
            ])

            ds.obs_sources = "\n".join([
                "http://opendap.aviso.oceanobs.com/",
                "http://hadobs.metoffice.com/",
                "http://oceandata.sci.gsfc.nasa.gov/"
            ])

            # ===== 변수 속성 (flag_values, flag_meanings) 추가 =====
            if 'obs_type' in ds.variables:
                var = ds.variables['obs_type']
                var.setncattr('flag_values', np.arange(1, 12, dtype='i4'))
                var.setncattr('flag_meanings',
                              "zeta ubar vbar u v temperature salinity NO3 phytoplankton zooplankton detritus")

            if 'obs_provenance' in ds.variables:
                var = ds.variables['obs_provenance']
                var.setncattr('flag_values', np.arange(1, 15, dtype='i4'))
                var.setncattr('flag_meanings',
                              "AVISO_SLA blended_SST XBT_temperature CTD_temperature CTD_salinity "
                              "ARGO_temperature ARGO_salinity CalCOFI_temperature CalCOFI_salinity "
                              "GLOBEC_temperature GLOBEC_salinity buoy_temp CTD_Chlorophyll Satellite_Chlorophyll")

        print(f"✅ Done: {ncfile}\n")


add_global_attrs_to_obs(outputFile)




