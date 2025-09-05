# -*- coding: utf-8 -*-
"""
Created on Thu May  8 18:00:19 2025

@author: ust21
"""
PKG_path = 'C:/Users/ust21/shjo/projects/myROMS/prc_src/' # Location of JNUROMS directory
import sys 
sys.path.append(PKG_path)
import utils.ROMS_utils01 as ru
import utils.ROMS_utils02 as ru2
from utils.ncCreate import create_ini_NPZD,create_std_NPZD2
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from tqdm import tqdm
from scipy.interpolate import griddata
from netCDF4 import Dataset,date2num,num2date

#== Define Inputs files =======================================================
My_Std='D:/shjo/ROMS_inputs/std/ROMS_std_15km_N20.nc' # Initial file name (to create)
My_Grd='D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc' # Grd name

#-- Define OGCM path ----------------------------------------------------------
ncdir='D:/shjo/ROMS_inputs/std/monthly_N20/'

#-- Define Parameters ---------------------------------------------------------
Ini_title='NWP 15km ROMS-NPZD std file' # title for NC description

# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
         'lon_u':'longitude','lat_u':'latitude','lon_v':'longitude','lat_v':'latitude',
         'temp':'temp_std','salt':'salt_std','u':'u_std','v':'v_std','zeta':'zeta_std',\
             'ubar':'ubar_std','vbar':'vbar_std'}

#== Starts Calc ===============================================================

#-- Get My Grid info ----------------------------------------------------------
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
MyVar={'Layer_N':20,'Vtransform':2,\
       'Vstretching':4,'Theta_s':6.5,\
           'Theta_b':1,'Tcline':400,'hmin':10}
# ncG.close()

create_std_NPZD2(ncdir+'ROMS_std_JAN_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_FEB_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_MAR_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_APR_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_MAY_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_JUN_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_JUL_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_AUG_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_SEP_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_OCT_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_NOV_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)
create_std_NPZD2(ncdir+'ROMS_std_DEC_15km_N20.nc',mask,topo,ncG,MyVar,0,Ini_title)


JAN=Dataset(ncdir+'ROMS_std_JAN_15km_N20.nc',mode='a')
FEB=Dataset(ncdir+'ROMS_std_FEB_15km_N20.nc',mode='a')
MAR=Dataset(ncdir+'ROMS_std_MAR_15km_N20.nc',mode='a')
APR=Dataset(ncdir+'ROMS_std_APR_15km_N20.nc',mode='a')
MAY=Dataset(ncdir+'ROMS_std_MAY_15km_N20.nc',mode='a')
JUN=Dataset(ncdir+'ROMS_std_JUN_15km_N20.nc',mode='a')
JUL=Dataset(ncdir+'ROMS_std_JUL_15km_N20.nc',mode='a')
AUG=Dataset(ncdir+'ROMS_std_AUG_15km_N20.nc',mode='a')
SEP=Dataset(ncdir+'ROMS_std_SEP_15km_N20.nc',mode='a')
OCT=Dataset(ncdir+'ROMS_std_OCT_15km_N20.nc',mode='a')
NOV=Dataset(ncdir+'ROMS_std_NOV_15km_N20.nc',mode='a')
DEC=Dataset(ncdir+'ROMS_std_DEC_15km_N20.nc',mode='a')

NC=Dataset('D:/shjo/ROMS_inputs/std/ROMS_std_15km_N20.nc')

# -----------------JAN---------------------------------------------------------
JAN['zeta'][0,:,:]=NC['zeta'][0,:,:]
JAN['ubar'][0,:,:]=NC['ubar'][0,:,:]
JAN['vbar'][0,:,:]=NC['vbar'][0,:,:]
JAN['temp'][0,:,:]=NC['temp'][0,:,:,:]
JAN['salt'][0,:,:]=NC['salt'][0,:,:,:]
JAN['u'][0,:,:]   =NC['u'][0,:,:,:]
JAN['v'][0,:,:]   =NC['v'][0,:,:,:] 
JAN['NO3'][0,:,:] =NC['NO3'][0,:,:,:]
JAN['phytoplankton'][0,:,:]=NC['phytoplankton'][0,:,:,:]
JAN['zooplankton'][0,:,:]  =NC['zooplankton'][0,:,:,:]
JAN['detritus'][0,:,:]     =NC['detritus'][0,:,:,:] 
JAN.close()

# -----------------FEB---------------------------------------------------------
FEB['zeta'][0,:,:]=NC['zeta'][1,:,:]
FEB['ubar'][0,:,:]=NC['ubar'][1,:,:]
FEB['vbar'][0,:,:]=NC['vbar'][1,:,:]
FEB['temp'][0,:,:]=NC['temp'][1,:,:,:]
FEB['salt'][0,:,:]=NC['salt'][1,:,:,:]
FEB['u'][0,:,:]   =NC['u'][1,:,:,:]
FEB['v'][0,:,:]   =NC['v'][1,:,:,:] 
FEB['NO3'][0,:,:] =NC['NO3'][1,:,:,:]
FEB['phytoplankton'][0,:,:]=NC['phytoplankton'][1,:,:,:]
FEB['zooplankton'][0,:,:]  =NC['zooplankton'][1,:,:,:]
FEB['detritus'][0,:,:]     =NC['detritus'][1,:,:,:] 
FEB.close()

# -----------------MAR---------------------------------------------------------
MAR['zeta'][0,:,:]=NC['zeta'][2,:,:]
MAR['ubar'][0,:,:]=NC['ubar'][2,:,:]
MAR['vbar'][0,:,:]=NC['vbar'][2,:,:]
MAR['temp'][0,:,:]=NC['temp'][2,:,:,:]
MAR['salt'][0,:,:]=NC['salt'][2,:,:,:]
MAR['u'][0,:,:]   =NC['u'][2,:,:,:]
MAR['v'][0,:,:]   =NC['v'][2,:,:,:] 
MAR['NO3'][0,:,:] =NC['NO3'][2,:,:,:]
MAR['phytoplankton'][0,:,:]=NC['phytoplankton'][2,:,:,:]
MAR['zooplankton'][0,:,:]  =NC['zooplankton'][2,:,:,:]
MAR['detritus'][0,:,:]     =NC['detritus'][2,:,:,:] 
MAR.close()

# -----------------APR---------------------------------------------------------
APR['zeta'][0,:,:]=NC['zeta'][3,:,:]
APR['ubar'][0,:,:]=NC['ubar'][3,:,:]
APR['vbar'][0,:,:]=NC['vbar'][3,:,:]
APR['temp'][0,:,:]=NC['temp'][3,:,:,:]
APR['salt'][0,:,:]=NC['salt'][3,:,:,:]
APR['u'][0,:,:]   =NC['u'][3,:,:,:]
APR['v'][0,:,:]   =NC['v'][3,:,:,:] 
APR['NO3'][0,:,:] =NC['NO3'][3,:,:,:]
APR['phytoplankton'][0,:,:]=NC['phytoplankton'][3,:,:,:]
APR['zooplankton'][0,:,:]  =NC['zooplankton'][3,:,:,:]
APR['detritus'][0,:,:]     =NC['detritus'][3,:,:,:] 
APR.close()

# -----------------MAY---------------------------------------------------------
MAY['zeta'][0,:,:]=NC['zeta'][4,:,:]
MAY['ubar'][0,:,:]=NC['ubar'][4,:,:]
MAY['vbar'][0,:,:]=NC['vbar'][4,:,:]
MAY['temp'][0,:,:]=NC['temp'][4,:,:,:]
MAY['salt'][0,:,:]=NC['salt'][4,:,:,:]
MAY['u'][0,:,:]   =NC['u'][4,:,:,:]
MAY['v'][0,:,:]   =NC['v'][4,:,:,:] 
MAY['NO3'][0,:,:] =NC['NO3'][4,:,:,:]
MAY['phytoplankton'][0,:,:]=NC['phytoplankton'][4,:,:,:]
MAY['zooplankton'][0,:,:]  =NC['zooplankton'][4,:,:,:]
MAY['detritus'][0,:,:]     =NC['detritus'][4,:,:,:] 
MAY.close()

# -----------------JUN---------------------------------------------------------
JUN['zeta'][0,:,:]=NC['zeta'][5,:,:]
JUN['ubar'][0,:,:]=NC['ubar'][5,:,:]
JUN['vbar'][0,:,:]=NC['vbar'][5,:,:]
JUN['temp'][0,:,:]=NC['temp'][5,:,:,:]
JUN['salt'][0,:,:]=NC['salt'][5,:,:,:]
JUN['u'][0,:,:]   =NC['u'][5,:,:,:]
JUN['v'][0,:,:]   =NC['v'][5,:,:,:] 
JUN['NO3'][0,:,:] =NC['NO3'][5,:,:,:]
JUN['phytoplankton'][0,:,:]=NC['phytoplankton'][5,:,:,:]
JUN['zooplankton'][0,:,:]  =NC['zooplankton'][5,:,:,:]
JUN['detritus'][0,:,:]     =NC['detritus'][5,:,:,:] 
JUN.close()

# -----------------JUL---------------------------------------------------------
JUL['zeta'][0,:,:]=NC['zeta'][6,:,:]
JUL['ubar'][0,:,:]=NC['ubar'][6,:,:]
JUL['vbar'][0,:,:]=NC['vbar'][6,:,:]
JUL['temp'][0,:,:]=NC['temp'][6,:,:,:]
JUL['salt'][0,:,:]=NC['salt'][6,:,:,:]
JUL['u'][0,:,:]   =NC['u'][6,:,:,:]
JUL['v'][0,:,:]   =NC['v'][6,:,:,:] 
JUL['NO3'][0,:,:] =NC['NO3'][6,:,:,:]
JUL['phytoplankton'][0,:,:]=NC['phytoplankton'][6,:,:,:]
JUL['zooplankton'][0,:,:]  =NC['zooplankton'][6,:,:,:]
JUL['detritus'][0,:,:]     =NC['detritus'][6,:,:,:] 
JUL.close()

# -----------------AUG---------------------------------------------------------
AUG['zeta'][0,:,:]=NC['zeta'][7,:,:]
AUG['ubar'][0,:,:]=NC['ubar'][7,:,:]
AUG['vbar'][0,:,:]=NC['vbar'][7,:,:]
AUG['temp'][0,:,:]=NC['temp'][7,:,:,:]
AUG['salt'][0,:,:]=NC['salt'][7,:,:,:]
AUG['u'][0,:,:]   =NC['u'][7,:,:,:]
AUG['v'][0,:,:]   =NC['v'][7,:,:,:] 
AUG['NO3'][0,:,:] =NC['NO3'][7,:,:,:]
AUG['phytoplankton'][0,:,:]=NC['phytoplankton'][7,:,:,:]
AUG['zooplankton'][0,:,:]  =NC['zooplankton'][7,:,:,:]
AUG['detritus'][0,:,:]     =NC['detritus'][7,:,:,:] 
AUG.close()

# -----------------SEP---------------------------------------------------------
SEP['zeta'][0,:,:]=NC['zeta'][8,:,:]
SEP['ubar'][0,:,:]=NC['ubar'][8,:,:]
SEP['vbar'][0,:,:]=NC['vbar'][8,:,:]
SEP['temp'][0,:,:]=NC['temp'][8,:,:,:]
SEP['salt'][0,:,:]=NC['salt'][8,:,:,:]
SEP['u'][0,:,:]   =NC['u'][8,:,:,:]
SEP['v'][0,:,:]   =NC['v'][8,:,:,:] 
SEP['NO3'][0,:,:] =NC['NO3'][8,:,:,:]
SEP['phytoplankton'][0,:,:]=NC['phytoplankton'][8,:,:,:]
SEP['zooplankton'][0,:,:]  =NC['zooplankton'][8,:,:,:]
SEP['detritus'][0,:,:]     =NC['detritus'][8,:,:,:] 
SEP.close()

# -----------------OCT---------------------------------------------------------
OCT['zeta'][0,:,:]=NC['zeta'][9,:,:]
OCT['ubar'][0,:,:]=NC['ubar'][9,:,:]
OCT['vbar'][0,:,:]=NC['vbar'][9,:,:]
OCT['temp'][0,:,:]=NC['temp'][9,:,:,:]
OCT['salt'][0,:,:]=NC['salt'][9,:,:,:]
OCT['u'][0,:,:]   =NC['u'][9,:,:,:]
OCT['v'][0,:,:]   =NC['v'][9,:,:,:] 
OCT['NO3'][0,:,:] =NC['NO3'][9,:,:,:]
OCT['phytoplankton'][0,:,:]=NC['phytoplankton'][9,:,:,:]
OCT['zooplankton'][0,:,:]  =NC['zooplankton'][9,:,:,:]
OCT['detritus'][0,:,:]     =NC['detritus'][9,:,:,:] 
OCT.close()

# -----------------NOV---------------------------------------------------------
NOV['zeta'][0,:,:]=NC['zeta'][10,:,:]
NOV['ubar'][0,:,:]=NC['ubar'][10,:,:]
NOV['vbar'][0,:,:]=NC['vbar'][10,:,:]
NOV['temp'][0,:,:]=NC['temp'][10,:,:,:]
NOV['salt'][0,:,:]=NC['salt'][10,:,:,:]
NOV['u'][0,:,:]   =NC['u'][10,:,:,:]
NOV['v'][0,:,:]   =NC['v'][10,:,:,:] 
NOV['NO3'][0,:,:] =NC['NO3'][10,:,:,:]
NOV['phytoplankton'][0,:,:]=NC['phytoplankton'][10,:,:,:]
NOV['zooplankton'][0,:,:]  =NC['zooplankton'][10,:,:,:]
NOV['detritus'][0,:,:]     =NC['detritus'][10,:,:,:] 
NOV.close()

# -----------------DEC---------------------------------------------------------
DEC['zeta'][0,:,:]=NC['zeta'][11,:,:]
DEC['ubar'][0,:,:]=NC['ubar'][11,:,:]
DEC['vbar'][0,:,:]=NC['vbar'][11,:,:]
DEC['temp'][0,:,:]=NC['temp'][11,:,:,:]
DEC['salt'][0,:,:]=NC['salt'][11,:,:,:]
DEC['u'][0,:,:]   =NC['u'][11,:,:,:]
DEC['v'][0,:,:]   =NC['v'][11,:,:,:] 
DEC['NO3'][0,:,:] =NC['NO3'][11,:,:,:]
DEC['phytoplankton'][0,:,:]=NC['phytoplankton'][11,:,:,:]
DEC['zooplankton'][0,:,:]  =NC['zooplankton'][11,:,:,:]
DEC['detritus'][0,:,:]     =NC['detritus'][11,:,:,:] 
DEC.close()














