####################################################################################################
# Author: Youngmin Park                                                                            #
# Description: This script reads ROMS result files and extracts variables from the standard area.  #
#              The extracted variables are saved to a new NetCDF file.                             #
# Date: 2025-03-20                                                                                 #
####################################################################################################


import numpy as np
import netCDF4 as nc
import yaml, os
import warnings
import sys
sys.path.append(os.path.join(os.getcwd(), 'libs'))
from libs.tools import getlist, ExtractVar
from libs.vgrid import compute_sigma_3d
from libs.interp_edit import LinearND_Nan, interpolate_parallel, vertical_interpolation_parallel
from libs.gennc_edit import GenStandardArea, create_netcdf_file, save_variable_to_netcdf, mkdir
import datetime, re

warnings.filterwarnings("ignore")


# open yaml file
with open('config.yml') as file:
    config = yaml.full_load(file)
    
    
RDomcfgPath = config['RDomcfgPath']
RomsRstHead = config['RomsRstHead']
RomsRstTail = config['RomsRstTail']
Romsnamlist = config['RomsNamlist']

output_head = config['OutputHead']

## config Header
cfgHeaders = config['nc_global_att']


rstlist = getlist(RomsRstHead, RomsRstTail)
for rstfile in rstlist:
    print("Read file: ", rstfile)
    
    
extvar = ExtractVar(Romsnamlist)
Vtransform = extvar.getvar('Vtransform')
Vstretching = extvar.getvar('Vstretching')
theta_s = extvar.getvar('THETA_S')
theta_b = extvar.getvar('THETA_B')
tcline  = extvar.getvar('TCLINE')
nlayer  = extvar.getvar('N')
## nlayer  = 30 print('!!! nlayer = 30 -> 36 !!!')
nlayer  = 36

# print vars
print('Vtransform  :', Vtransform)
print('Vstretching :', Vstretching)
print('theta_s     :', theta_s)
print('theta_b     :', theta_b)
print('tcline      :', tcline)
print('nlayer      :', nlayer)


romsdomcfg = nc.Dataset(RDomcfgPath, 'r')
print("Read Domain file: ", RDomcfgPath)


romlat_rho  = romsdomcfg.variables['lat_rho'][:]
romlon_rho  = romsdomcfg.variables['lon_rho'][:]
romlat_u    = romsdomcfg.variables['lat_u'][:]
romlon_u    = romsdomcfg.variables['lon_u'][:]
romlat_v    = romsdomcfg.variables['lat_v'][:]
romlon_v    = romsdomcfg.variables['lon_v'][:]
romlat_psi  = romsdomcfg.variables['lat_psi'][:]
romlon_psi  = romsdomcfg.variables['lon_psi'][:]

rommask_rho = romsdomcfg.variables['mask_rho'][:]
rommask_u   = romsdomcfg.variables['mask_u'][:]
rommask_v   = romsdomcfg.variables['mask_v'][:]
rommask_psi = romsdomcfg.variables['mask_psi'][:]
romh        = romsdomcfg.variables['h'][:]


print('lat_rho shape:', romlat_rho.shape)
print('lon_rho shape:', romlon_rho.shape)
print('lat_u   shape:', romlat_u.shape)
print('lon_u   shape:', romlon_u.shape)
print('lat_v   shape:', romlat_v.shape)
print('lon_v   shape:', romlon_v.shape)
print('lat_psi shape:', romlat_psi.shape)
print('lon_psi shape:', romlon_psi.shape)



stdlatmin = config['Standard']['latmin']
stdlatmax = config['Standard']['latmax']
stdlonmin = config['Standard']['lonmin']
stdlonmax = config['Standard']['lonmax']
stdhresol = config['Standard']['hresol']
stdvdepth = config['Standard']['vdepth']
savehead  = config['StdSaveHead']

mkdir(savehead)

print('Standard Area:')
print('latmin:', stdlatmin)
print('latmax:', stdlatmax)
print('lonmin:', stdlonmin)
print('lonmax:', stdlonmax)
print('hresol:', stdhresol)
print('vdepth:', stdvdepth)


std_area, std_lat_grid, std_lon_grid, std_depth = GenStandardArea(stdlatmin, stdlatmax, stdlonmin, stdlonmax, stdhresol, stdvdepth)


print('Standard Area shape:', std_area.shape)
print('Standard lat grid shape:', std_lat_grid.shape)
print('Standard lon grid shape:', std_lon_grid.shape)
print('Standard depth shape:', std_depth.shape)

## 기존 코드
# std_h   = LinearND_Nan(romlon_rho, romlat_rho, romh, std_lon_grid, std_lat_grid)
# depth3d = compute_sigma_3d(std_h, tcline, theta_s, theta_b, nlayer, Vtransform, Vstretching) * -1
# depth3d = depth3d[::-1, :, :]
# print('depth3d shape:', depth3d.shape)

stdvars = config['ExtraVar']



for rstpath in rstlist:
    
    rst = nc.Dataset(rstpath, 'r')
    
    ## 기준 시간 단위, 시작 시간
    strlist = [x.strip() for x in rst.variables['ocean_time'].units.split('since')]
    ref_timeunit = strlist[0]
    ref_date_str = strlist[1]
    
    ## 시간 계산
    target_ocean_time = rst.variables['ocean_time'][:][0]
    target_ocean_time = (target_ocean_time * 60 if 'minute' in ref_timeunit else 
                         target_ocean_time * 3600 if 'hour' in ref_timeunit else 
                         target_ocean_time * 86400 if 'day' in ref_timeunit else 
                         target_ocean_time)
    date_obj = datetime.datetime.strptime(ref_date_str, '%Y-%m-%d %H:%M:%S')
    date_obj = date_obj + datetime.timedelta(seconds=target_ocean_time) # - datetime.timedelta(days=1)
    date_str = date_obj.strftime('%Y%m%d')
    
    ## 폴더명 반영
    # stdncpath = os.path.join(savehead, 'std_' + os.path.basename(rstpath))
    # stdncpath = re.sub(r'(\D+)_\d+(\.nc)', r'\1_' + date_str + r'\2', stdncpath)
    stdncpath = savehead + output_head + date_str + '.nc' 
    print('\nProcessing start :', stdncpath)
    
    stdnc = create_netcdf_file(stdncpath, std_lat_grid, std_lon_grid, std_depth, cfgHeaders)
    
    
    ## Zeta 반영
    zeta = rst.variables['zeta'][:]  # (time, lat, lon) 형태로 가정
    zeta = np.squeeze(zeta)  # 시간 차원 제거 (필요 시)
    
    std_h   = LinearND_Nan(romlon_rho, romlat_rho, romh + zeta, std_lon_grid, std_lat_grid)
    depth3d = compute_sigma_3d(std_h, tcline, theta_s, theta_b, nlayer, Vtransform, Vstretching) * -1
    depth3d = depth3d[::-1, :, :]
    

    for var_dict in stdvars:
        var_name = list(var_dict.keys())[0]
        grid_type, vartype = var_dict[var_name]
        
        print(f' - Processing variable: {var_name} (Grid: {grid_type}, Type: {vartype})')

        roms_var = np.squeeze(rst.variables[var_name][:])
        if len(roms_var.shape) ==3:
            roms_var = roms_var[::-1, :, :]

        if grid_type == 'u':
            romlon, romlat = romlon_u, romlat_u
        elif grid_type == 'v':
            romlon, romlat = romlon_v, romlat_v
        else:
            romlon, romlat = romlon_rho, romlat_rho

        roms_var_interpolated = interpolate_parallel(romlon, romlat, roms_var, std_lon_grid, std_lat_grid, n_jobs=-1)

        if vartype == 3:
            std_array = vertical_interpolation_parallel(depth3d, std_h, roms_var_interpolated, std_depth, n_jobs=-1)
        else:
            std_array = roms_var_interpolated  # Directly store 2D data
        
        save_variable_to_netcdf(stdnc, var_name, std_array, vartype)
        
        
    
    print("Done processing: ", stdncpath)

    stdnc.close()
    rst.close()
