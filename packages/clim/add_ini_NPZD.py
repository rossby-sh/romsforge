# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 13:04:04 2025

@author: ust21
"""

PKG_path = '/home/shjo/github/romsforge/libs/' # Location of JNUROMS directory
import sys 
sys.path.append(PKG_path)
import ROMS_utils01 as ru
import ROMS_utils02 as ru2
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from tqdm import tqdm
from scipy.interpolate import griddata
from netCDF4 import Dataset,date2num,num2date

#== Define Inputs files =======================================================
My_Ini='/home/shjo/applications/nifs02_clm/data/eccov4_clm5km_ini_fennel.nc' # Initial file name (to create)
My_Grd='/home/shjo/data/roms_inputs/grd/mcc/roms_grd_fennel_5km_smooth_v3.nc' # Grd name
#-- Define OGCM path ----------------------------------------------------------

ncdir='/home/shjo/data/raw/cmems_bio_hist/clm_CMEMS_bio.nc'
NO3NC=ncdir
phytNC=ncdir
#o2NC=ncdir+'BIO/CMEMS_data_bio_2025-02.nc'




#-- Define Parameters ---------------------------------------------------------
Ini_title='my Test ROMS' # title for NC description

conserv=1
# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
         'NO3':'no3','phyt':'phyc','zoop':'zooc','detr':'???'}

# Define time info

#== Starts Calc ===============================================================

#-- Get My Grid info ----------------------------------------------------------
ncG=Dataset(My_Grd)

print(ncG)

lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
MyVar={'Layer_N':36,'Vtransform':2,\
       'Vstretching':4,'Theta_s':6.5,\
           'Theta_b':1,'Tcline':400,'hmin':10}
ncG.close()

atG,onG=lonG.shape
cosa,sina=np.cos(angle),np.sin(angle)

#-- Get OGCM Grid info --------------------------------------------------------
ncO=Dataset(NO3NC)
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

#-- Create a dump ncfile ------------------------------------------------------
TIMES_co=0
#-- Get OGCM data for initial -------------------------------------------------
OGCM_Data={}#,OGCM_Mask={}
for i in ['NO3','phyt']:
#for i in ['temp']:

    print('!!! Data processing : '+i+' !!!')

    if i=='NO3' :
        OGCM_npth=NO3NC;
    elif i=='phyt':
        OGCM_npth=phytNC;
    elif i=='zoop':
        OGCM_npth=zoopNC;
    elif i=='detr':
        OGCM_npth=detrNC;
   
    tmp_data=np.squeeze(Dataset(OGCM_npth)[OGCMVar[i]][TIMES_co,:,latO_co,lonO_co])
    tmp_mask=np.invert(tmp_data.mask)
## Another way to create mask
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
    
OGCM_Data['detr']=np.ones_like(OGCM_Data['phyt'])*0.04
OGCM_Data['zoop']=OGCM_Data['phyt']*0.3

    
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


ncI=Dataset(My_Ini,mode='a')
ncI['NO3'][0]=NO3
ncI['phytoplankton'][0]=phyt
ncI['zooplankton'][0]=zoop
ncI['detritus'][0]=detr
ncI.close()







