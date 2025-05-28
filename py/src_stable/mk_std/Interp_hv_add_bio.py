# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 13:04:04 2025

@author: ust21
"""

PKG_path = 'C:/Users/ust21/shjo/projects/myROMS/prc_src/' # Location of JNUROMS directory
import sys 
sys.path.append(PKG_path)
import utils.ROMS_utils01 as ru
import utils.ROMS_utils02 as ru2
from utils.ncCreate import create_ini_NPZD
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
STDNC=ncdir+'glorys_21y_Bio_std.nc'

#-- Define Parameters ---------------------------------------------------------

# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
         'NO3':'no3_std','phyt':'phyc_std','zoop_std':'zooc','detr':'???'}

#-- Get My Grid info ----------------------------------------------------------
ncG=Dataset(My_Grd)

print(ncG)

lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
MyVar={'Layer_N':20,'Vtransform':2,\
       'Vstretching':4,'Theta_s':6.5,\
           'Theta_b':1,'Tcline':400,'hmin':10}
ncG.close()

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

#-- Create a dump ncfile ------------------------------------------------------
# create_ini_NPZD(My_Ini,mask,topo,MyVar,Ini_time_num,Ini_title,ncFormat='NETCDF4')

#-- Get OGCM data for initial -------------------------------------------------
OGCM_Data={}#,OGCM_Mask={}

for mm in range(12):
    for i in ['NO3','phyt']:
    #for i in ['temp']:
    
        print('!!! Data processing : '+i+' !!!')

       
        tmp_data=np.squeeze(Dataset(STDNC)[OGCMVar[i]][mm,:,latO_co,lonO_co])
        #tmp_mask=np.invert(tmp_data.mask)
    ## Another way to create mask
        tmp_mask=tmp_data!=tmp_data
        #tmp_mask=np.invert(tmp_data.mask)
        tmp_mask=np.invert(tmp_mask)
    
    
    # mv=ncO[OGCMVar[i]].missing_value
    
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
        
    OGCM_Data['detr']=OGCM_Data['phyt']*0.1
    OGCM_Data['zoop']=OGCM_Data['phyt']*0.1

        
    #-- Process vector elements ---------------------------------------------------
    
    #-- Process ROMS Vertical grid ------------------------------------------------
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000
    
    _,y_tmp,x_tmp=OGCM_Data['NO3'].shape
    Rzeta=np.zeros([y_tmp,x_tmp])
    
    zr=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 1,topo,Rzeta);
        
    zw=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 5,topo,Rzeta);
    dzr=zw[1:,:,:]-zw[:-1,:,:]
    
    
    #-- Add a level on top and bottom with no-gradient ----------------------------
    NO3,phyt=OGCM_Data['NO3'],OGCM_Data['phyt']
    zoop,detr=OGCM_Data['zoop'],OGCM_Data['detr']
    
    NO3=np.vstack((np.expand_dims(NO3[0,:,:],axis=0)\
                  ,NO3,np.expand_dims(NO3[-1,:,:],axis=0)))
    phyt=np.vstack((np.expand_dims(phyt[0,:,:],axis=0)\
                  ,phyt,np.expand_dims(phyt[-1,:,:],axis=0)))
    zoop=np.vstack((np.expand_dims(zoop[0,:,:],axis=0)\
                  ,zoop,np.expand_dims(zoop[-1,:,:],axis=0)))
    detr=np.vstack((np.expand_dims(detr[0,:,:],axis=0)\
                  ,detr,np.expand_dims(detr[-1,:,:],axis=0)))
    
    #-- Transform z-coordinate to sigma-coordinates -------------------------------
    print('!!! Transformming z --> sigma... !!!')
    
    NO3=ru.ztosigma(np.flip(NO3,axis=0),zr,np.flipud(Z));
    phyt=ru.ztosigma(np.flip(phyt,axis=0),zr,np.flipud(Z));
    zoop=ru.ztosigma(np.flip(zoop,axis=0),zr,np.flipud(Z));
    detr=ru.ztosigma(np.flip(detr,axis=0),zr,np.flipud(Z));
    #-- Volume conservation -------------------------------------------------------
        
    #-- Calc Barotropic velocities2 -----------------------------------------------
    
    ncI=Dataset(My_Std,mode='a')
    ncI['NO3'][mm,:,:,:]=NO3
    ncI['phytoplankton'][mm,:,:,:]=phyt
    ncI['zooplankton'][mm,:,:,:]=zoop
    ncI['detritus'][mm,:,:,:]=detr
    ncI.close()
    
    
