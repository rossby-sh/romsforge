# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 09:44:00 2024

@author: shjo
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
ncdir='D:/shjo/ROMS_inputs/std/'
STDNC=ncdir+'glorys_21y_Phy_std.nc'


#-- Define Parameters ---------------------------------------------------------
Ini_title='NWP 15km ROMS-NPZD std file' # title for NC description

# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
         'lon_u':'longitude','lat_u':'latitude','lon_v':'longitude','lat_v':'latitude',
         'temp':'temp_std','salt':'salt_std','u':'u_std','v':'v_std','zeta':'zeta_std',\
             'ubar':'ubar_std','vbar':'vbar_std'}

# Define time info
t_rng=['2025-01-01','2025-01-01'] # Inital time 
My_time_ref='days since 2000-1-1 00:00:00' # time ref

#== Starts Calc ===============================================================

#-- Get My Grid info ----------------------------------------------------------
ncG=Dataset(My_Grd)

print(ncG)

lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
MyVar={'Layer_N':20,'Vtransform':2,\
       'Vstretching':4,'Theta_s':6.5,\
           'Theta_b':1,'Tcline':400,'hmin':10}
# ncG.close()

atG,onG=lonG.shape
cosa,sina=np.cos(angle),np.sin(angle)

#-- Get OGCM Grid info --------------------------------------------------------
ncO=Dataset(STDNC)
lonO,latO=ncO[OGCMVar['lon_rho']][:],ncO[OGCMVar['lat_rho']][:];
depthO=ncO[OGCMVar['depth']][:]
ncO.close()

#-- Get OGCM lon lat coordinates for slicing ----------------------------------

lonO_co01=np.where( (lonO[:]>=np.min(lonG)) & (lonO[:]<=np.max(lonG)) )[0]
latO_co01=np.where( (latO[:]>=np.min(latG)) & (latO[:]<=np.max(latG)) )[0]

latO_re=latO[latO_co01]
lonO_re=lonO[lonO_co01]

lonGmax=np.max(np.abs(np.diff(lonG[0,:])))
latGmax=np.max(np.abs(np.diff(latG[:,0])))

lonOmax=np.max(np.abs(np.diff(lonO_re[:])))
latOmax=np.max(np.abs(np.diff(latO_re[:])))

lonEval=np.max([ np.max(lonO_re) ,np.max(lonG)+lonOmax])
lonSval=np.min([ np.min(lonO_re), np.min(lonG)-lonOmax])

latEval=np.max([ np.max(latO_re),np.max(latG)+latOmax])
latSval=np.min([ np.min(latO_re),np.min(latG)-latOmax])

lonO_co=np.where((lonO[:]>=lonSval)&(lonO[:]<=lonEval))[0]
latO_co=np.where((latO[:]>=latSval)&(latO[:]<=latEval))[0]

latO_s=latO[latO_co]
lonO_s=lonO[lonO_co]

lonO_s_m,latO_s_m=np.meshgrid(lonO_s,latO_s)


#-- Process Times--------------------------------------------------------------
# OGCM_TIMES=Dataset(STDNC)[OGCMVar['time']]
# TIME_UNIT=OGCM_TIMES.units
# OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
# Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),int(t_rng[0].split('-')[2]),0)
# Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),int(t_rng[1].split('-')[2]),23)
# TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# tmp_y,tmp_m,tmp_d=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[1]),int(t_rng[1].split('-')[-1])
# tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,tmp_d),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,tmp_d),My_time_ref)
# Ini_time_num=((date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-tmp_dif))
# # print('!!! Ini_time + 16d !!!')

# #-- Create a dump ncfile ------------------------------------------------------

Ini_time_num=np.arange(15.5,366,30.4375)*86400

create_std_NPZD2(My_Std,mask,topo,ncG,MyVar,Ini_time_num,Ini_title)

#-- Get OGCM data for initial -------------------------------------------------
OGCM_Data={}#,OGCM_Mask={}

for mm in range(12):
        
    for i in ['zeta','ubar','vbar','u','v','temp','salt']:
    # for i in ['zeta']:
    
        print('!!! Data processing : '+i+' !!!')
        if (i == 'zeta') or (i == 'ubar') or (i == 'vbar'):
            tmp_data=np.squeeze(Dataset(STDNC)[OGCMVar[i]][mm,latO_co,lonO_co])
            tmp_mask=tmp_data!=tmp_data
            # tmp_mask=np.invert(tmp_data.mask)
            tmp_mask=np.invert(tmp_mask)
        
        else:
           
            tmp_data=np.squeeze(Dataset(STDNC)[OGCMVar[i]][mm,:,latO_co,lonO_co])
            
            tmp_mask=tmp_data!=tmp_data
            #tmp_mask=np.invert(tmp_data.mask)
            tmp_mask=np.invert(tmp_mask)
            
        ## Another way to create mask
        # mv=ncO[OGCMVar[i]].missing_value
        if len(tmp_data.shape)==2:
    
            data=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                          tmp_data[tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
            data_ex=data.reshape(latO_s_m.shape)
            data=griddata((lonO_s_m.flatten(),latO_s_m.flatten()),\
                          data_ex.flatten(),(lonG.flatten(),latG.flatten()),'cubic')
            data=data.reshape(lonG.shape)
            tmp_var=data
            if np.sum(np.isnan(data))!=0:
                tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                data=tmp_var 
           
        elif len(tmp_data.shape)==3:
        
            data,n=np.zeros([len(depthO),atG,onG]),0
            for j,k in tqdm(zip(tmp_data,tmp_mask)):
         
                data_=griddata((lonO_s_m[k].flatten(),latO_s_m[k].flatten()),\
                              j[k].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_ex=data_.reshape(latO_s_m.shape)
                data_=griddata((lonO_s_m.flatten(),latO_s_m.flatten()),\
                              data_ex.flatten(),(lonG.flatten(),latG.flatten()),'cubic')
                data[n]=data_.reshape(latG.shape)
                
            ### Add 20250320
                tmp_var=data[n]
                if np.sum(~np.isnan(data[n]))==0:
                    data[n]=data[n-1]
    
                if np.sum(np.isnan(data[n]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[n]=tmp_var
                n+=1
        OGCM_Data[i]=data
        
    #-- Process vector elements ---------------------------------------------------
    
    ubar_=OGCM_Data['ubar']
    vbar_=OGCM_Data['vbar']
    
    ubar=ru2.rho2u_2d(ubar_)
    vbar=ru2.rho2v_2d(vbar_)

    
    #-- Process ROMS Vertical grid ------------------------------------------------
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000
    
    Rzeta=np.zeros_like(OGCM_Data['zeta'])
    
    
    zr=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 1,topo,Rzeta);
    zu=ru2.rho2u_3d(zr);
    zv=ru2.rho2v_3d(zr);
    zw=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 5,topo,Rzeta);
    dzr=zw[1:,:,:]-zw[:-1,:,:]
    dzu=ru2.rho2u_3d(dzr);
    dzv=ru2.rho2v_3d(dzr);
    
    #-- Add a level on top and bottom with no-gradient ----------------------------
    temp,salt=OGCM_Data['temp'],OGCM_Data['salt']
    u,v=OGCM_Data['u'],OGCM_Data['v']
    
    u1=np.vstack((np.expand_dims(u[0,:,:],axis=0)\
                  ,u,np.expand_dims(u[-1,:,:],axis=0)))
    v1=np.vstack((np.expand_dims(v[0,:,:],axis=0)\
                  ,v,np.expand_dims(v[-1,:,:],axis=0)))
    temp=np.vstack((np.expand_dims(temp[0,:,:],axis=0)\
                  ,temp,np.expand_dims(temp[-1,:,:],axis=0)))
    salt=np.vstack((np.expand_dims(salt[0,:,:],axis=0)\
                  ,salt,np.expand_dims(salt[-1,:,:],axis=0)))
    
    #-- Transform z-coordinate to sigma-coordinates -------------------------------
    print('!!! Transformming z --> sigma... !!!')
    u=ru.ztosigma(np.flip(u1,axis=0),zu,np.flipud(Z));
    v=ru.ztosigma(np.flip(v1,axis=0),zv,np.flipud(Z));
    temp=ru.ztosigma(np.flip(temp,axis=0),zr,np.flipud(Z));
    salt=ru.ztosigma(np.flip(salt,axis=0),zr,np.flipud(Z));
    
    #-- Calc Barotropic velocities2 -----------------------------------------------
    ncI=Dataset(My_Std,mode='a')
    ncI['zeta'][mm,:,:]=OGCM_Data['zeta']
    ncI['ubar'][mm,:,:]=ubar
    ncI['vbar'][mm,:,:]=vbar
    # ncI['SSH'][:]=Rzeta
    ncI['temp'][mm,:,:,:]=temp
    ncI['salt'][mm,:,:,:]=salt
    ncI['u'][mm,:,:,:]=u
    ncI['v'][mm,:,:,:]=v
    
    ncI.close()







