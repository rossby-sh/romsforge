import numpy as np
import netCDF4 as nc

def GenStandardArea(latmin, latmax, lonmin, lonmax, hresol, vdepth):

    lat = np.linspace(latmin, latmax, int((latmax - latmin) / hresol) + 1)
    lon = np.linspace(lonmin, lonmax, int((lonmax - lonmin) / hresol) + 1)

    lon_grid, lat_grid = np.meshgrid(lon, lat)

    depth = np.array(vdepth)

    Nz = len(depth)
    Ny, Nx = lat_grid.shape
    area_3d = np.zeros((Nz, Ny, Nx))

    return area_3d, lat_grid, lon_grid, depth




def create_netcdf_file(output_path, std_lat_grid, std_lon_grid, std_depth):

    lat_size, lon_size = std_lat_grid.shape
    depth_size = len(std_depth)

    ncfile = nc.Dataset(output_path, 'w', format='NETCDF4')

    ncfile.createDimension('lat', lat_size)
    ncfile.createDimension('lon', lon_size)
    ncfile.createDimension('depth', depth_size)
    ncfile.createDimension('time', None)

    nc_lat = ncfile.createVariable('lat', 'f4', ('lat', 'lon'))
    nc_lon = ncfile.createVariable('lon', 'f4', ('lat', 'lon'))
    nc_depth = ncfile.createVariable('depth', 'f4', ('depth',))

    nc_lat[:, :] = std_lat_grid
    nc_lon[:, :] = std_lon_grid
    nc_depth[:] = std_depth

    nc_lat.units = "degrees_north"
    nc_lon.units = "degrees_east"
    nc_depth.units = "meters"
    nc_depth.positive = "down"

    return ncfile




def save_variable_to_netcdf(ncfile, var_name, var_data, vartype):

    if vartype == 2:
        if var_name not in ncfile.variables:
            ncfile.createVariable(var_name, 'f4', ('lat', 'lon'))
        ncfile.variables[var_name][:] = var_data

    elif vartype == 3:
        if var_name not in ncfile.variables:
            ncfile.createVariable(var_name, 'f4', ('depth', 'lat', 'lon'))
        ncfile.variables[var_name][:] = var_data

    ncfile.sync()



def mkdir(path):
    import os
    if not os.path.exists(path):
        print(f"Creating directory: {path}")
        os.makedirs(path)