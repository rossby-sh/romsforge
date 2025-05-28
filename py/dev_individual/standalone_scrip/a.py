from netCDF4 import Dataset

f = Dataset('remap_grid_fennel_15km_smooth_v2_rho.nc')

print("✅ Variables:", list(f.variables.keys()))
print("✅ Dimensions:", dict(f.dimensions))
print("✅ grid_dims dtype:", f.variables['grid_dims'].dtype)
print("✅ grid_center_lon dtype:", f.variables['grid_center_lon'].dtype)

