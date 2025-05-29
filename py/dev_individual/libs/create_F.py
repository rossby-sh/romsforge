# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:25:43 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np


def createF_era5(nc_save_name, LON, LAT, TIME, Ref_time, variables, ncFormat='NETCDF3_CLASSIC'):
    """
    Create a ROMS-style forcing NetCDF file using ERA5-like variables.

    Parameters:
    - nc_save_name: output NetCDF file name
    - LON, LAT: 1D longitude and latitude arrays
    - TIME: 1D time array (shared by all variables)
    - Ref_time: time reference string (e.g., 'days since 2000-01-01')
    - variables: dict of {varname: 3D ndarray (time, lat, lon)}
    """
 
    lat_len = len(LAT)
    lon_len = len(LON)

    forcing_variables = {
        'Uwind':       ('f4', ('wind_time', 'lat', 'lon'), {'long_name': 'Eastward Wind', 'units': 'm/s'}),
        'Vwind':       ('f4', ('wind_time', 'lat', 'lon'), {'long_name': 'Northward Wind', 'units': 'm/s'}),
        'Tair':        ('f4', ('tair_time', 'lat', 'lon'), {'long_name': 'Air Temperature', 'units': 'Celsius'}),
        'Qair':        ('f4', ('qair_time', 'lat', 'lon'), {'long_name': 'Specific Humidity', 'units': 'kg/kg'}),
        'Cloud':       ('f4', ('cloud_time', 'lat', 'lon'), {'long_name': 'Cloud Fraction', 'units': 'percent'}),
        'sst':         ('f4', ('sst_time', 'lat', 'lon'), {'long_name': 'Sea Surface Temperature', 'units': 'Celsius'}),
        'dqdsst':      ('f4', ('dqdsst_time', 'lat', 'lon'), {'long_name': 'dQ/dSST', 'units': 'Watts/m^2/C'}),
        'srf':         ('f4', ('srf_time', 'lat', 'lon'), {'long_name': 'Net Shortwave Radiation Flux', 'units': 'Watts/m^2'}),
        'lwrad':       ('f4', ('lwrad_time', 'lat', 'lon'), {'long_name': 'Net Longwave Radiation Flux', 'units': 'Watts/m^2'}),
        'lwrad_down':  ('f4', ('lrad_down_time', 'lat', 'lon'), {'long_name': 'Downward Longwave Radiation', 'units': 'Watts/m^2'}),
        'rain':        ('f4', ('rain_time', 'lat', 'lon'), {'long_name': 'Precipitation Rate', 'units': 'kg/m^2/s'}),
        'Pair':        ('f4', ('pair_time', 'lat', 'lon'), {'long_name': 'Surface Air Pressure', 'units': 'millibars'}),
    }

    time_metadata = {
        time_name: {
            'long_name': f'{time_name.replace("_time", "").capitalize()} forcing time',
            'units': Ref_time
        }
        for _, dims, _ in forcing_variables.values()
        for time_name in dims if time_name.endswith('_time')
    }

    ncfile = Dataset(nc_save_name, mode='w', format=ncFormat)

    # Define spatial dimensions
    ncfile.createDimension('lon', lon_len)
    ncfile.createDimension('lat', lat_len)
    for time_name in time_metadata:
        dim_len = None if ncFormat == 'NETCDF4' else len(TIME)
        ncfile.createDimension(time_name, dim_len)

    # Coordinate variables
    lon = ncfile.createVariable('lon', 'f4', ('lon',))
    lat = ncfile.createVariable('lat', 'f4', ('lat',))
    lon.long_name = 'Longitude'
    lon.units = 'degrees_east'
    lat.long_name = 'Latitude'
    lat.units = 'degrees_north'

    # Time variables (define only)
    time_vars = {}
    for time_name, attrs in time_metadata.items():
        tvar = ncfile.createVariable(time_name, 'f8', (time_name,))
        tvar.long_name = attrs['long_name']
        tvar.units = attrs['units']
        tvar.field = f"{time_name}, scalar, series"
        time_vars[time_name] = tvar

    var_objects = {}
    # Forcing variables
    for varname, data in variables.items():
        if varname not in forcing_variables:
            raise KeyError(f"Unknown variable: {varname}")

        dtype, dims, attrs = forcing_variables[varname]
        v = ncfile.createVariable(varname, dtype, dims)
        for attr, val in attrs.items():
            setattr(v, attr, val)
        v.coordinates = 'lon lat'
        v.time = dims[0]  # 첫 번째 차원이 시간
        v.field = f"{varname}, scalar, series"
        var_objects[varname] = (v, data)
    
    # Now assign time values after all variables are defined
    for tname in time_vars:
        time_vars[tname][:] = TIME

    for varname, (v, data) in var_objects.items():
        v[:,:,:] = data

    lon[:] = LON
    lat[:] = LAT

    ncfile.title = 'ROMS ERA5-style Forcing'
    ncfile.history = 'Created using create_forcing_era5'

    ncfile.close()
    print(f"✅ NetCDF saved: {nc_save_name}")



def createF_era5_(nc_save_name, LON, LAT, TIME, Ref_time, variables,ncFormat='NETCDF3_CLASSIC'):
    """
    Create a single NetCDF file with multiple forcing variables.

    Parameters:
    - nc_save_name: output NetCDF file name
    - LON, LAT: 1D longitude and latitude arrays
    - TIME: 1D time array
    - Ref_time: time reference string, e.g. 'days since 2000-01-01'
    - variables: dictionary where keys are variable names (e.g., 'Uwind') and values are 3D numpy arrays (time, lat, lon)
    """
    ncfile = Dataset(nc_save_name, mode='w', format=ncFormat)

    # Create dimensions
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('time', len(TIME))

    # Create coordinate variables
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    time = ncfile.createVariable('time', np.float64, ('time',))

    lon.long_name = 'Longitude'
    lon.units = 'degrees_east'
    lat.long_name = 'Latitude'
    lat.units = 'degrees_north'
    time.long_name = 'Time'
    time.units = Ref_time
    time.field = 'time, scalar, series'

    lon[:] = LON
    lat[:] = LAT
    time[:] = TIME

    # Define variable metadata (can be customized per variable)
    metadata = {
        'Uwind': dict(long_name='Eastward Wind', units='m/s'),
        'Vwind': dict(long_name='Northward Wind', units='m/s'),
        'Tair': dict(long_name='Air Temperature', units='Celsius'),
        'Qair': dict(long_name='Specific Humidity', units='kg/kg'),
        'Cloud': dict(long_name='Cloud Fraction', units='percent'),
        'sst': dict(long_name='Sea Surface Temperature', units='Celsius'),
        'dqdsst': dict(long_name='dQ/dSST', units='Watts/m^2/C'),
        'srf': dict(long_name='Net Shortwave Radiation Flux', units='Watts/m^2'),
        'lwrad': dict(long_name='Net Longwave Radiation Flux', units='Watts/m^2'),
        'lwrad_down': dict(long_name='Downward Longwave Radiation', units='Watts/m^2'),
        'rain': dict(long_name='Precipitation Rate', units='kg/m^2/s'),
        'Pair': dict(long_name='Surface Air Pressure', units='millibars'),
    }

    # Write all variables
    for varname, data in variables.items():
        var = ncfile.createVariable(varname, np.float64, ('time', 'lat', 'lon'))
        var[:, :, :] = data

        # Apply standard metadata
        var.coordinates = 'lon lat'
        var.time = 'time'
        if varname in metadata:
            var.long_name = metadata[varname].get('long_name', varname)
            var.units = metadata[varname].get('units', 'unknown')
        else:
            var.long_name = varname
            var.units = 'unknown'

    ncfile.close()
    print(f"✅ NetCDF saved: {nc_save_name}")

def get_dqdsst(sst,sat,rho_atm,U,qsea):
    '''
    %  sst     : sea surface temperature (Celsius)
    %  sat     : sea surface atmospheric temperature (Celsius)
    %  rho_atm : atmospheric density (kilogram meter-3) 
    %  U       : wind speed (meter s-1)
    %  qsea    : sea level specific humidity
    '''
    # Specific heat of atmosphere.
    Cp=1004.8
    # Sensible heat transfert coefficient (stable condition)
    Ch = 0.66e-3
    # Latent heat transfert coefficient (stable condition)
    Ce = 1.15e-3
    # Emissivity coefficient
    eps = 0.98
    # Stefan constant
    stef = 5.6697e-8;
    # SST (KELVIN)
    SST = sst + 273.15;
    # Latent heat of vaporisation (J.kg-1)
    L = 2.5008e6 - 2.3e3 * sat
    # Infrared contribution
    q1 = -4. * stef * (SST**3)
    # Sensible heat contribution
    q2 = -rho_atm * Cp * Ch * U
    # Latent heat contribution
    dqsdt = 2353.* np.log(10.) * qsea / (SST**2)
    q3 = -rho_atm * Ce * L * U * dqsdt
    dqdsst = q1 + q2 + q3 
    return dqdsst

def create_wind_nc(nc_save_name,LON,LAT,TIME,Ref_time,values1,values2):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('wind_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly winds'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('wind_time', np.float64, ('wind_time',))
    time.units=Ref_time
    time.field='qair_time, scalar, series'
    # time.cycle_length=cycle
    
    Uwind = ncfile.createVariable('Uwind',np.float64,('wind_time','lat','lon')) # note: unlimited dimension is leftmost
    Uwind.units = 'meter second-1' 
    Uwind.long_name = 'surface u-wind component' # this is a CF standard name
    Uwind.time='wind_time'
    Uwind.coordinates = "lon lat"
    Uwind.field='u-wind, scalar, series'
    
    Vwind = ncfile.createVariable('Vwind',np.float64,('wind_time','lat','lon')) # note: unlimited dimension is leftmost
    Vwind.units = 'meter second-1' 
    Vwind.long_name = 'surface v-wind component' # this is a CF standard name
    Vwind.time='wind_time'
    Vwind.coordinates = "lon lat"
    Vwind.field='v-wind, scalar, series'


    # Data.field=Var.field

    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Uwind[:] = values1
    Vwind[:] = values2

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
    
def create_tair_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('tair_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly tair'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('tair_time', np.float64, ('tair_time',))
    time.units=Ref_time
    time.field='tair_time, scalar, series'
    Tair = ncfile.createVariable('Tair',np.float64,('tair_time','lat','lon')) # note: unlimited dimension is leftmost
    Tair.units = 'Celsius' 
    Tair.long_name = 'surface air temperature at 2m, for CORE' # this is a CF standard name
    Tair.time='tair_time'
    Tair.coordinates = "lon lat"
    Tair.field='Tair, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Tair[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
        
def create_pair_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('pair_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly pair (1993)'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('pair_time', np.float64, ('pair_time',))
    time.units=Ref_time
    time.field='pair_time, scalar, series'

    
    Pair = ncfile.createVariable('Pair',np.float64,('pair_time','lat','lon')) # note: unlimited dimension is leftmost
    Pair.units = 'millibar' 
    Pair.long_name = 'surface air pressure, for CORE' # this is a CF standard name
    Pair.time='pair_time'
    Pair.coordinates = "lon lat"
    Pair.field='Tair, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Pair[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    

def create_qair_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('qair_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly qair'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('qair_time', np.float64, ('qair_time',))
    time.units=Ref_time
    # time.cycle_length=cycle
    time.field='qair_time, scalar, series'

    
    Qair = ncfile.createVariable('Qair',np.float64,('qair_time','lat','lon')) # note: unlimited dimension is leftmost
    Qair.units = 'percentage' 
    Qair.long_name = 'surface air specific humidity, for CORE' # this is a CF standard name
    Qair.time='qair_time'
    Qair.coordinates = "lon lat"
    Qair.field='Qair, scalar, series'
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Qair[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
    
def create_cloud_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('cloud_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly cloud'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('cloud_time', np.float64, ('cloud_time',))
    time.units=Ref_time
    # time.cycle_length=cycle
    time.field='cloud_time, scalar, series'

    
    Cloud = ncfile.createVariable('Cloud',np.float64,('cloud_time','lat','lon')) # note: unlimited dimension is leftmost
    Cloud.units = 'nondimensional' 
    Cloud.long_name = '-' # this is a CF standard name
    Cloud.time='cloud_time'
    Cloud.coordinates = "lon lat"
    Cloud.field='Cloud, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Cloud[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    

def create_sst_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('sst_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly sst'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('sst_time', np.float64, ('sst_time',))
    time.units=Ref_time
    time.field='sst_time, scalar, series'

    
    SST = ncfile.createVariable('sst',np.float64,('sst_time','lat','lon')) # note: unlimited dimension is leftmost
    SST.units =  'Celsius'
    SST.long_name = 'sea surface temperature' # this is a CF standard name
    SST.time='sst_time'
    SST.coordinates = "lon lat"
    SST.field='SST, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    SST[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
def create_dqdsst_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('sst_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly dqdsst '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('sst_time', np.float64, ('sst_time',))
    time.units=Ref_time
    time.field='sst_time, scalar, series'

    
    dqdsst = ncfile.createVariable('dqdsst',np.float64,('sst_time','lat','lon')) # note: unlimited dimension is leftmost
    dqdsst.units = 'Watts meter-2 Celsius-1' 
    dqdsst.long_name = 'surface net heat flux sensitivity to SST' # this is a CF standard name
    dqdsst.time='sst_time'
    dqdsst.coordinates = "lon lat"
    dqdsst.field='dQdSST, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    dqdsst[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
def create_srf_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('srf_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly srf '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('srf_time', np.float64, ('srf_time',))
    time.units=Ref_time
    time.field='srf_time, scalar, series'

    
    srf = ncfile.createVariable('srf',np.float64,('srf_time','lat','lon')) # note: unlimited dimension is leftmost
    srf.units = 'Watt meter-2'
    srf.long_name = 'shortwave radiation, scalar, series' # this is a CF standard name
    srf.time='srf_time'
    srf.positive='downward flux, heating'
    srf.negative='upward flux, cooling'
    srf.coordinates = "lon lat"
    srf.field='shortwave radiation, scalar, series'


    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    srf[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
def create_lwrad_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('lrf_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly lwrad'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('lrf_time', np.float64, ('lrf_time',))
    time.units=Ref_time
    time.field='lrf_time, scalar, series'

    
    lrf = ncfile.createVariable('lwrad',np.float64,('lrf_time','lat','lon')) # note: unlimited dimension is leftmost
    lrf.units = 'Watt meter-2'
    lrf.long_name = 'Net longwave radiation' # this is a CF standard name
    lrf.time='lrf_time'
    lrf.positive='downward flux, heating'
    lrf.negative='upward flux, cooling'
    lrf.coordinates = "lon lat"
    lrf.field = 'longwave radiation, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    lrf[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()   
    
    
def create_lwrad_down_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('lrf_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly lwrad_down'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('lrf_time', np.float64, ('lrf_time',))
    time.units=Ref_time
    time.field='lrf_time, scalar, series'

    lwrad_down = ncfile.createVariable('lwrad_down',np.float64,('lrf_time','lat','lon')) # note: unlimited dimension is leftmost
    lwrad_down.units = 'Watt meter-2'
    lwrad_down.long_name = 'Downward longwave radiation' # this is a CF standard name
    lwrad_down.time='lrf_time'
    lwrad_down.positive='downward flux, heating'
    lwrad_down.negative='upward flux, cooling'
    lwrad_down.coordinates = "lon lat"
    lwrad_down.field = 'longwave radiation, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    lwrad_down[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()   
    
    
    
   
def create_rain_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('rain_time',len(TIME))
    
    ncfile.title='ERA5 3 hourly rain'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('rain_time', np.float64, ('rain_time',))
    time.units=Ref_time
    time.field='rain_time, scalar, series'

    
    rain = ncfile.createVariable('rain',np.float64,('rain_time','lat','lon')) # note: unlimited dimension is leftmost
    rain.units = 'kg mm-2 s-1'
    rain.long_name = 'rain fall rate' # this is a CF standard name
    rain.time='rain_time'
    rain.coordinates = "lon lat"
    rain.field = 'rain, scalar, series'

    
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    rain[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()   
    
    
    
    
    
    
