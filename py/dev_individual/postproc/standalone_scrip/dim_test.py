from netCDF4 import Dataset

file_paths = [
    "remap_grid_GLBv0.08_Arctic4_t.nc",
    "remap_grid_fennel_15km_smooth_v2_rho.nc"
]

for path in file_paths:
    with Dataset(path) as nc:
        dims = nc.variables.get("grid_dims")
        print("\nFile:", path)
        if dims is not None:
            print("  dtype:", dims.dtype)
            print("  value:", dims[:])
        else:
            print("  grid_dims not found")

