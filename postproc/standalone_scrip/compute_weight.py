from pyroms.remapping import compute_remap_weights

compute_remap_weights(
    grid1_file='remap_grid_GLBv0.08_Arctic4_t.nc',  # HYCOM (source)
    grid2_file='remap_grid_fennel_15km_smooth_v2_rho.nc',  # ROMS (destination)
    interp_file1='weights_HYCOM_to_ROMS_bilinear_temp.nc',
    interp_file2='weights_ROMS_to_HYCOM_bilinear_temp.nc',
    map1_name='HYCOM to ROMS Bilinear Mapping',
    map2_name='ROMS to HYCOM Bilinear Mapping',
    num_maps=2,
    map_method='bilinear'
)
