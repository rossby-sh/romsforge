
import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d
from joblib import Parallel, delayed
from concurrent.futures import ThreadPoolExecutor


def LinearND_Nan(lon1, lat1, data, lon2, lat2) -> np.ndarray:
    
    if len(np.shape(lon1)) == 1 or len(np.shape(lat1)) == 1:
        lon1, lat1 = np.meshgrid(lon1, lat1)
    
    if len(np.shape(lon2)) == 1 or len(np.shape(lat2)) == 1:
        lon2, lat2 = np.meshgrid(lon2, lat2)

    data = np.where(data.mask, np.nan, data)

    Ndinterp = LinearNDInterpolator((lon1.ravel(), 
                                     lat1.ravel()), 
                                     data.ravel())

    filled_val = Ndinterp(lon2, lat2)

    return filled_val



def interpolate_parallel(romlon_rho, romlat_rho, roms_var, std_lon_grid, std_lat_grid, n_jobs=-1):

    if len(roms_var.shape) == 2:
        horibox = LinearND_Nan(romlon_rho, romlat_rho, roms_var, std_lon_grid, std_lat_grid)
        return horibox
    
    elif len(roms_var.shape) == 3:
        horibox = np.zeros((roms_var.shape[0], std_lon_grid.shape[0], std_lon_grid.shape[1]))

        results = Parallel(n_jobs=n_jobs)(
            delayed(LinearND_Nan)(romlon_rho, romlat_rho, roms_var[i, :, :], std_lon_grid, std_lat_grid) 
            for i in range(roms_var.shape[0])
        )
        horibox[:] = np.array(results)
        return horibox

    elif len(roms_var.shape) == 4:
        raise ValueError("4D input not supported yet")

    else:
        raise ValueError("Input data must be 2D, 3D, or 4D")
    
    

def vertical_interpolation_parallel(depth3d, roms_var, stdvdepth, n_jobs=-1):
    
    num_std_layers = len(stdvdepth)
    lat_size, lon_size = depth3d.shape[1], depth3d.shape[2]
    std_area = np.full((num_std_layers, lat_size, lon_size), np.nan)

    def interpolate_column(i):
        interpolated_column = np.full((num_std_layers, lat_size), np.nan)
        for j in range(lat_size):
            depth_profile = depth3d[:, j, i]
            value_profile = roms_var[:, j, i]

            valid_mask = ~np.isnan(depth_profile) & ~np.isnan(value_profile)
            if np.sum(valid_mask) > 2:
                interp_func = interp1d(depth_profile[valid_mask], value_profile[valid_mask], 
                                       kind='linear', bounds_error=False, 
                                       fill_value=(value_profile[valid_mask][0], value_profile[valid_mask][-1]))  
                interpolated_column[:, j] = interp_func(stdvdepth)
        return interpolated_column

    results = Parallel(n_jobs=n_jobs)(delayed(interpolate_column)(i) for i in range(lon_size))

    std_area[:] = np.stack(results, axis=-1)

    return std_area