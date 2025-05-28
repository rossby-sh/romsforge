# remap_module.py

from netCDF4 import Dataset
import numpy as np
from scipy.sparse import coo_matrix
import os
import xarray as xr
import xesmf as xe
import glob


def read_grid(grid_nc):
    ds = xr.open_dataset(grid_nc)
    coord_vars = {}
    for name in ['lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v']:
        if name in ds:
            coord_vars[name] = ds[name]
    return coord_vars


def create_regridder(src_grid, dst_grid, method="bilinear", weight_path=None):
    regridder = xe.Regridder(src_grid, dst_grid, method=method, filename=weight_path, reuse_weights=False)
    regridder.to_netcdf(weight_path)
    return regridder


def load_weight_matrix(weight_nc):
    nc = Dataset(weight_nc)
    row = nc.variables['row'][:]
    col = nc.variables['col'][:]
    S = nc.variables['S'][:]
    dst_size = row.max() + 1
    src_size = col.max() + 1
    W = coo_matrix((S, (row, col)), shape=(dst_size, src_size))
    return W


def apply_weight_remap_2d(weight_matrix, src_2d, dst_shape):
    src_flat = src_2d.ravel()
    dst_flat = weight_matrix.dot(src_flat)
    return dst_flat.reshape(dst_shape)


def apply_weight_remap_3d(weight_matrix, src_3d, dst_shape):
    Nz = src_3d.shape[0]
    dst_3d = np.empty((Nz, *dst_shape), dtype=src_3d.dtype)
    for k in range(Nz):
        dst_3d[k] = apply_weight_remap_2d(weight_matrix, src_3d[k], dst_shape)
    return dst_3d


def duplicate_ncStruct(input_file, output_file, new_coords):
    src = Dataset(input_file)
    dst = Dataset(output_file, 'w')

    for name, dim in src.dimensions.items():
        dst.createDimension(name, (len(dim) if not dim.isunlimited() else None))

    for name, var in src.variables.items():
        if name in new_coords:
            continue
        dst_var = dst.createVariable(name, var.datatype, var.dimensions, zlib=True)
        for attr in var.ncattrs():
            dst_var.setncattr(attr, var.getncattr(attr))
        dst_var[:] = 0

    for name, arr in new_coords.items():
        if arr.ndim == 2:
            dims = ('eta_' + name[-1], 'xi_' + name[-1])
        else:
            continue
        dst.createVariable(name, 'f4', dims)[:] = arr

    dst.close()
    src.close()


def process_files(avg_dir, grid01, grid02, weight_paths, output_dir):
    print("[INFO] Reading grids...")
    src_coords = read_grid(grid01)
    dst_coords = read_grid(grid02)

    grids = {}
    weights = {}

    for key in ['rho', 'u', 'v']:
        lon_key = f'lon_{key}'
        lat_key = f'lat_{key}'
        if lon_key in src_coords and lat_key in src_coords:
            src_grid = xr.Dataset({'lat': src_coords[lat_key], 'lon': src_coords[lon_key]})
            dst_grid = xr.Dataset({'lat': dst_coords[lat_key], 'lon': dst_coords[lon_key]})
            dst_shape = dst_coords[lat_key].shape
            weight_nc = weight_paths[key]

            if not os.path.exists(weight_nc):
                print(f"[INFO] Creating weight file for {key}-points...")
                create_regridder(src_grid, dst_grid, "bilinear", weight_nc)

            grids[key] = dst_shape
            weights[key] = load_weight_matrix(weight_nc)

    avg_files = sorted(glob.glob(os.path.join(avg_dir, '*.nc')))
    for file in avg_files:
        print(f"[INFO] Processing {file}...")
        base = os.path.basename(file)
        out_file = os.path.join(output_dir, base)

        duplicate_ncStruct(file, out_file, dst_coords)

        src = Dataset(file)
        dst = Dataset(out_file, 'a')

        for var in src.variables:
            dims = src.variables[var].dimensions
            if 'xi_rho' in dims:
                key = 'rho'
            elif 'xi_u' in dims:
                key = 'u'
            elif 'xi_v' in dims:
                key = 'v'
            else:
                continue

            data = src.variables[var][:]
            dst_shape = grids[key]
            W = weights[key]

            if data.ndim == 3:
                remapped = apply_weight_remap_3d(W, data, dst_shape)
            elif data.ndim == 2:
                remapped = apply_weight_remap_2d(W, data, dst_shape)
            else:
                continue
            dst.variables[var][:] = remapped

        dst.close()
        src.close()
        print(f"[INFO] Saved remapped file: {out_file}")
