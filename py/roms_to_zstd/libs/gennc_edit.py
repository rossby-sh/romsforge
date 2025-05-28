import numpy as np
import netCDF4 as nc
import datetime

def GenStandardArea(latmin, latmax, lonmin, lonmax, hresol, vdepth):

    lat = np.linspace(latmin, latmax, int((latmax - latmin) / hresol) + 1)
    lon = np.linspace(lonmin, lonmax, int((lonmax - lonmin) / hresol) + 1)

    lon_grid, lat_grid = np.meshgrid(lon, lat)

    depth = np.array(vdepth)

    Nz = len(depth)
    Ny, Nx = lat_grid.shape
    area_3d = np.zeros((Nz, Ny, Nx))

    return area_3d, lat_grid, lon_grid, depth




def create_netcdf_file(output_path, std_lat_grid, std_lon_grid, std_depth, headers):

    lat_size, lon_size = std_lat_grid.shape
    depth_size = len(std_depth)
    
    ncformat='NETCDF4'

    ncfile = nc.Dataset(output_path, 'w', format=ncformat)

    ncfile.createDimension('lat', lat_size)
    ncfile.createDimension('lon', lon_size)
    ncfile.createDimension('depth', depth_size)
    #ncfile.createDimension('time', None)

    nc_lat = ncfile.createVariable('lat', 'f4', ('lat', 'lon'))
    nc_lon = ncfile.createVariable('lon', 'f4', ('lat', 'lon'))
    nc_depth = ncfile.createVariable('depth', 'f4', ('depth',))

    nc_lat[:, :] = std_lat_grid
    nc_lon[:, :] = std_lon_grid
    nc_depth[:] = std_depth

    nc_lon.standard_name = "grid_longitude_at_cell_center" ;
    nc_lon.long_name = "longitude of RHO-points" ;
    nc_lon.units = "degrees_east"

    nc_lat.standard_name = "grid_latitude_at_cell_center" ;
    nc_lat.long_name = "latitude of RHO-points" ;
    nc_lat.units = "degrees_north"

    nc_depth.standard_name = "standard_depth" ;
    nc_depth.long_name = "standard_depth" ;
    nc_depth.units = "meters"
    nc_depth.positive = "down"
    
    ncfile.format = ncformat
    ncfile.type = headers['type']
    ncfile.model_time_ref = headers['model_time_ref']
    ncfile.created_by = headers['created_by']
    ncfile.history = datetime.datetime.now().strftime('%Y.%m.%d %H:%M') + ' (KST)'

    return ncfile




def save_variable_to_netcdf(ncfile, var_name, var_data, vartype):

    if vartype == 2:
        if var_name not in ncfile.variables:
            nc_tmp = ncfile.createVariable(var_name, 'f4', ('lat', 'lon'))
            
            if var_name == 'zeta':
                nc_tmp.standard_name = "sea_surface_height_above_geopotential_datum"
                nc_tmp.long_name = "free-surface"
                nc_tmp.units = "meter"
                
        ncfile.variables[var_name][:] = var_data

    elif vartype == 3:
        if var_name not in ncfile.variables:
            nc_tmp = ncfile.createVariable(var_name, 'f4', ('depth', 'lat', 'lon'))
            if var_name == 'temp':
                nc_tmp.standard_name = "sea_water_potential_temperature"
                nc_tmp.long_name = "potential temperature"
                nc_tmp.units = "Celsius"
                
            elif var_name == 'salt':
                nc_tmp.standard_name = "sea_water_practical_salinity"
                nc_tmp.long_name = "salinity"
                
            elif var_name == 'u':
                nc_tmp.standard_name = "sea_water_x_velocity"
                nc_tmp.long_name = "u-momentum component"
                nc_tmp.units = "meter second-1"
                
            elif var_name == 'v':
                nc_tmp.standard_name = "sea_water_y_velocity"
                nc_tmp.long_name = "v-momentum component"
                nc_tmp.units = "meter second-1"
            
            elif var_name == 'NO3':
                nc_tmp.standard_name = "mole_concentration_of_nitrate_expressed_as_nitrogen_in_sea_water"
                nc_tmp.long_name = "nitrate concentration"
                nc_tmp.units = "millimole_nitrogen meter-3"

            elif var_name == 'chlorophyll':
                nc_tmp.standard_name = "mole_concentration_of_chlorophyll_in_sea_water"
                nc_tmp.long_name = "chlorophyll concentration"
                nc_tmp.units = "milligrams_chlorophyll meter-3"
                
            elif var_name == 'phytoplankton':
                nc_tmp.standard_name = "mole_concentration_of_phytoplankton_expressed_as_nitrogen_in_sea_water"
                nc_tmp.long_name = "phytoplankton concentration"
                nc_tmp.units = "millimole_nitrogen meter-3"
                
            elif var_name == 'zooplankton':
                nc_tmp.standard_name = "mole_concentration_of_zooplankton_expressed_as_nitrogen_in_sea_water"
                nc_tmp.long_name = "zooplankton concentration"
                nc_tmp.units = "millimole_nitrogen meter-3"

        ncfile.variables[var_name][:] = var_data

    ncfile.sync()



def mkdir(path):
    import os
    if not os.path.exists(path):
        print(f"Creating directory: {path}")
        os.makedirs(path)