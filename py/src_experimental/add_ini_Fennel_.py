# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 13:04:04 2025

@author: ust21
"""

PKG_path = '/home/shjo/ROMS/romsforge/py/dev_individual/' # Location of JNUROMS directory
import sys 
sys.path.append(PKG_path)
import libs.ROMS_utils01 as ru
import libs.ROMS_utils02 as ru2
import numpy as np
import datetime as dt
from tqdm import tqdm
from scipy.interpolate import griddata
from netCDF4 import Dataset,date2num,num2date

#== Define Inputs files =======================================================
My_Ini='/data/share/DATA/ROMS_INPUTS/tmp/NWP4_ini_3_10m_LP.nc' # Initial file name (to create)
My_Grd='/data/share/DATA/ROMS_INPUTS/grd/NWP4_grd_3_10m_LP.nc' # Grd name

#-- Define OGCM path ----------------------------------------------------------
ncdir='/data/share/DATA/RAW/Bvar/'
NO3NC=ncdir+'NUT/CMEMS_data_nut_2023-01.nc'
phytNC=ncdir+'PFT/CMEMS_data_pft_2023-01.nc'
o2NC=ncdir+'BIO/CMEMS_data_bio_2023-01.nc'

#-- Define Parameters ---------------------------------------------------------

# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
        'NO3':'no3','PO4':'po4','chlorophyll':'chl','oxygen':'o2','phytoplankton':'phyc'}

# Define time info
t_rng=['2023-01-01','2023-01-01'] # Inital time 
My_time_ref='seconds since 2000-1-1 00:00:00' # time ref

#== Starts Calc ===============================================================

#-- Get My Grid info ----------------------------------------------------------
ncG=Dataset(My_Grd)

print(ncG)

lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
MyVar={'Layer_N':20,'Vtransform':2,\
       'Vstretching':2,'Theta_s':7,\
           'Theta_b':0.1,'Tcline':200,'hmin':10}
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


#-- Process Times--------------------------------------------------------------
OGCM_TIMES=Dataset(NO3NC)[OGCMVar['time']]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),int(t_rng[0].split('-')[2]),0)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),int(t_rng[1].split('-')[2]),23)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
tmp_y,tmp_m,tmp_d=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[1]),int(t_rng[1].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,tmp_d),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,tmp_d),My_time_ref)
Ini_time_num=((date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-tmp_dif))*86400
# print('!!! Ini_time + 16d !!!')

# debugging
print('=== time_co ===')
print(TIMES_co)
#-- Create a dump ncfile ------------------------------------------------------
# create_ini_NPZD(My_Ini,mask,topo,MyVar,Ini_time_num,Ini_title,ncFormat='NETCDF4')

#-- Get OGCM data for initial -------------------------------------------------
OGCM_Data={}#,OGCM_Mask={}
for i in ['NO3','phytoplankton','PO4','chlorophyll','oxygen']:
#for i in ['temp']:

    print('!!! Data processing : '+i+' !!!')

    if (i=='NO3') or (i=='PO4'):
        OGCM_npth=NO3NC;
    elif (i=='phytoplankton') or (i=='chlorophyll'):
        OGCM_npth=phytNC;
    elif i=='oxygen':
        OGCM_npth=o2NC;
   
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
NO3,PO4=OGCM_Data['NO3'],OGCM_Data['PO4']
chl,phyt=OGCM_Data['chlorophyll'],OGCM_Data['phytoplankton']
oxygen=OGCM_Data['oxygen']


NO3=np.vstack((np.expand_dims(NO3[0,:,:],axis=0)\
              ,NO3,np.expand_dims(NO3[-1,:,:],axis=0)))
PO4=np.vstack((np.expand_dims(PO4[0,:,:],axis=0)\
              ,PO4,np.expand_dims(PO4[-1,:,:],axis=0)))
chl=np.vstack((np.expand_dims(chl[0,:,:],axis=0)\
              ,chl,np.expand_dims(chl[-1,:,:],axis=0)))
phyt=np.vstack((np.expand_dims(phyt[0,:,:],axis=0)\
              ,phyt,np.expand_dims(phyt[-1,:,:],axis=0)))
oxygen=np.vstack((np.expand_dims(oxygen[0,:,:],axis=0)\
              ,oxygen,np.expand_dims(oxygen[-1,:,:],axis=0)))

#-- Transform z-coordinate to sigma-coordinates -------------------------------
print('!!! Transformming z --> sigma... !!!')

NO3=ru.ztosigma(np.flip(NO3,axis=0),zr,np.flipud(Z));
PO4=ru.ztosigma(np.flip(PO4,axis=0),zr,np.flipud(Z));
chl=ru.ztosigma(np.flip(chl,axis=0),zr,np.flipud(Z));
phyt=ru.ztosigma(np.flip(phyt,axis=0),zr,np.flipud(Z));
oxygen=ru.ztosigma(np.flip(oxygen,axis=0),zr,np.flipud(Z));
#-- Volume conservation -------------------------------------------------------

    
#-- Calc Barotropic velocities2 -----------------------------------------------


ncI=Dataset(My_Ini,mode='a')
ncI['NO3'][0]=NO3
ncI['PO4'][0]=PO4
ncI['NH4'][0]=0
ncI['TIC'][0]=2100
ncI['oxygen'][0]=oxygen
ncI['chlorophyll'][0]=chl
ncI['phytoplankton'][0]=phyt
ncI['zooplankton'][0]=phyt*0.3
ncI['alkalinity'][0]=2350
ncI['SdetritusC'][0]=0
ncI['LdetritusC'][0]=0
ncI['RdetritusC'][0]=0
ncI['SdetritusN'][0]=0
ncI['LdetritusN'][0]=0
ncI['RdetritusN'][0]=0
ncI.close()
















