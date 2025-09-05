import pyroms
import pyroms_toolbox

# remapping 모듈이 scrip을 인식하고 있어야 함
import scrip
pyroms.remapping.scrip = scrip

from pyroms.remapping import compute_remap_weights

# 입력 파일 이름
grid1_file = 'remap_grid_GLBv0.08_Arctic4_t.nc'
grid2_file = 'remap_grid_fennel_15km_smooth_v2_rho.nc'
interp_file1 = 'weights_HYCOM_to_ROMS_bilinear_temp.nc'
interp_file2 = 'weights_ROMS_to_HYCOM_bilinear_temp.nc'
map1_name = 'HYCOM to ROMS Bilinear Mapping'
map2_name = 'ROMS to HYCOM Bilinear Mapping'

# remap 수행
compute_remap_weights(
    grid1_file=grid1_file,
        grid2_file=grid2_file,
            interp_file1=interp_file1,
                interp_file2=interp_file2,
                    map1_name=map1_name,
                        map2_name=map2_name,
                            num_maps=1,
                                map_method='bilinear'
                                )

