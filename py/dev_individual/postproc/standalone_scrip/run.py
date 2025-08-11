import pyroms
import pyroms_toolbox
from pyroms_toolbox.Grid_HYCOM import get_nc_Grid_HYCOM2

# HYCOM netCDF 경로
h_npth = '/home/shjo/ROMS/_data/HYCOM/HYCOM4NWP12_TOT_230103_03.nc4'

# ROMS grid 객체
dst_grd = pyroms.grid.get_ROMS_grid('fennel_15km_smooth_v2')  # gridid.txt 확인

# HYCOM grid 객체
src_grd = get_nc_Grid_HYCOM2(h_npth)

# remap grid 생성
pyroms_toolbox.Grid_HYCOM.make_remap_grid_file(src_grd)
pyroms.remapping.make_remap_grid_file(dst_grd, Cpos='rho')
