# -*- coding: utf-8 -*-
"""
Created on Mon May 19 13:25:15 2025

@author: ust21
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
import datetime as dt
My_time_ref='day since 2000-1-1 00:00:00' # time ref

# num2date(9224,My_time_ref)
# date2num(dt.datetime(2025,4,3),My_time_ref)

mynum=9224
npth='D:/shjo/ROMS_inputs/NWP15km/N36/ROMS_bry_15km_250228_250501_N36.nc'

BRY=Dataset(npth,'a')

TIME=BRY['bry_time'][:]

t_co=np.where(TIME==mynum)[0]

BRY.variables.keys()

my_keys=['u_north', 'v_north', 'temp_north', 'salt_north', 'ubar_north', 'vbar_north', 'zeta_north',\
     'phyt_north', 'NO3_north', 'zoop_north', 'detritus_north',\
         'u_south', 'v_south', 'temp_south', 'salt_south', 'ubar_south', 'vbar_south', 'zeta_south',\
             'phyt_south', 'NO3_south', 'zoop_south', 'detritus_south',\
                 'u_east', 'v_east', 'temp_east', 'salt_east', 'ubar_east', 'vbar_east', 'zeta_east',\
                     'phyt_east', 'NO3_east', 'zoop_east', 'detritus_east']

for nm in my_keys:
    
    if (nm.split('_')[0] == 'zeta') or (nm.split('_')[0] == 'ubar') or (nm.split('_')[0] == 'vbar'):
        print('!!! '+nm)
        BRY[nm][t_co,:]=BRY[nm][t_co-1,:]
    
    else:
        BRY[nm][t_co,:,:]=BRY[nm][t_co-1,:,:]
        print('### '+nm)
    
BRY.close()
    
    