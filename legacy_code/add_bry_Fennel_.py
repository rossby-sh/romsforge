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
import xarray as xr
import os 
My_Bry='/data/shjo/nifs01/NWP12_bry_edit.nc' # Initial file name (to create)
My_Grd='/data/shjo/nifs01/NWP12_grd_edit_depth.nc' # Grd name
 
Parallel=False

#== Define Inputs files =======================================================

#-- Define OGCM path ----------------------------------------------------------
#-- Define OGCM path ----------------------------------------------------------
ncdir='/data/share/DATA/RAW/Bvar/'
NO3NC=ncdir+'NUT/'
phytNC=ncdir+'PFT/'
o2NC=ncdir+'BIO/'


NSEW=[1,1,1,1]

# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
        'NO3':'no3','PO4':'po4','chlorophyll':'chl','oxygen':'o2','phytoplankton':'phyc'}

conserv=1

OGCMS=[ncdir+'NUT/'+i for i in os.listdir(ncdir+'NUT/') if (i.endswith('.nc') & i.startswith('CMEMS_') ) ]

# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]

MyVar={'Layer_N':36,'Vtransform':2,\
       'Vstretching':2,'Theta_s':7,\
           'Theta_b':0.1,'Tcline':200,'hmin':10}
ncG.close()


# Get OGCM Grid info
ncO=Dataset(OGCMS[0])

lonO,latO=ncO[OGCMVar['lon_rho']][:],ncO[OGCMVar['lat_rho']][:]
depthO=ncO[OGCMVar['depth']][:]

#t_rng=['2019-12-01 00:00','2019-12-02 23:00']
#My_time_ref='days since 1990-1-1 00:00:00'
OGCM_TIMES=xr.open_mfdataset(OGCMS,decode_times=False)[OGCMVar['time']]

#OGCM_TIMES=MFDataset(OGCMS)[OGCMVar['time']] 
TIME_UNIT=OGCM_TIMES.units

t_rng = ['2022-12-31 00:00', '2025-01-01 23:00']
# t_rng = ['2024-12-31 00:00', '2025-04-01 23:00']
My_time_ref = 'days since 1980-01-01 00:00:00'
TIME_UNIT = OGCM_TIMES.units

# Convert OGCM time to datetime (e.g. HYCOM time is in "seconds since 1970-01-01")
OGCM_times = num2date(OGCM_TIMES[:], TIME_UNIT)

# Define target datetime range
Tst = dt.datetime.strptime(t_rng[0], "%Y-%m-%d %H:%M")
Ted = dt.datetime.strptime(t_rng[1], "%Y-%m-%d %H:%M")

# Get the time indices in range
TIMES_co = np.where((OGCM_times >= Tst) & (OGCM_times <= Ted))[0]
print(TIMES_co)
# Convert to new reference time base
OGCM_days = date2num(OGCM_times[TIMES_co], TIME_UNIT)  # time in original unit
Bry_time_num = date2num(OGCM_times[TIMES_co], My_time_ref)  # re-referenced time
Bry_time_time = num2date(Bry_time_num, My_time_ref)

thO=ncO[OGCMVar['depth']].shape[0]
ncO.close()

atG,onG=lonG.shape


#-- Get OGCM lon lat coordinates for slicing ----------------------------------
lonO_co01=np.where( (lonO[:]>=np.min(lonG)) & (lonO[:]<=np.max(lonG)) )[0]
latO_co01=np.where( (latO[:]>=np.min(latG)) & (latO[:]<=np.max(latG)) )[0]

latO_re=latO[latO_co01]
lonO_re=lonO[lonO_co01]

lonGmax=np.max(np.abs(np.diff(lonG[:])))
latGmax=np.max(np.abs(np.diff(latG[:])))

lonOmax=np.max(np.abs(np.diff(lonO_re[:])))
latOmax=np.max(np.abs(np.diff(latO_re[:])))

lonEval=np.max([ np.max(lonO_re) ,np.max(lonG)+lonOmax+0.5])
lonSval=np.min([ np.min(lonO_re), np.min(lonG)-lonOmax-0.5])

latEval=np.max([ np.max(latO_re),np.max(latG)+latOmax+0.5])
latSval=np.min([ np.min(latO_re),np.min(latG)-latOmax-0.5])

lonO_co=np.where((lonO[:]>=lonSval)&(lonO[:]<=lonEval))[0]
latO_co=np.where((latO[:]>=latSval)&(latO[:]<=latEval))[0]

latO_s=latO[latO_co]
lonO_s=lonO[lonO_co]

# lonO_s_m,latO_s_m=np.meshgrid(lonO_s,latO_s)
# lonO_s_m,latO_s_m=np.meshgrid(lonO_s,latO_s)

# =============================================================================
# Northern Boundary
# =============================================================================
if NSEW[0]:
    bry_lat_co=np.where( (latO[:]>=np.min(latG[-2:,:])-latOmax) & (latO[:]<=np.max(latG[-2:,:])+latOmax) )[0]
    bry_lon_co=lonO_co
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

    OGCM_Data={}#,OGCM_Mask={}
    for i in ['NO3','phytoplankton','PO4','chlorophyll','oxygen']:
    #for i in ['temp']:

        print('!!! Data processing : '+i+' !!!')

        if (i=='NO3') or (i=='PO4'):
            OGCM_npth=NO3NC;
            print('!!1')
        elif (i=='phytoplankton') or (i=='chlorophyll'):
            OGCM_npth=phytNC;
            print('!!!2')
        elif i=='oxygen':
            OGCM_npth=o2NC;
            print('!!!3')

       # tmp_data=np.squeeze(MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co])
        print(OGCM_npth+'*.nc')
        tmp_data=xr.open_mfdataset(OGCM_npth+'*.nc',parallel=Parallel)[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
        #tmp_data=MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co]  
        
        mask = np.isnan(tmp_data)
        DATA = np.ma.array(tmp_data, mask=mask)
        
        tmp_mask_=np.invert(DATA.mask)    
        data=np.zeros([len(Bry_time_num),len(depthO),2,lonG.shape[-1]])
    
        for t in tqdm(range(len(Bry_time_num))):
            for d in range(len(depthO)):
                # Interp mask
                tmp_mask=tmp_mask_[t,d]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                                DATA[t,d][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
                
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[-2:,:],latG[-2:,:]) ,'cubic' )
                data[t,d]=data_re_.reshape(lonG[-2:,:].shape) #.reshape(-1)
                
                tmp_var=data[t,d]
                
                if np.sum(~np.isnan(data[t,d]))==0:
                    data[t,d]=data[t,d-1]

                if np.sum(np.isnan(data[t,d]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t,d]=tmp_var
                    
        OGCM_Data[i]=data
        
    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    _,_,y_tmp,x_tmp=OGCM_Data['NO3'].shape

    Rzeta=np.zeros([len(TIMES_co),y_tmp,x_tmp])


    zr_= np.zeros([len(TIMES_co),MyVar['Layer_N'],Rzeta.shape[1],Rzeta.shape[-1]])
    zw = np.zeros([len(TIMES_co),MyVar['Layer_N']+1,Rzeta.shape[1],Rzeta.shape[-1]])

    
    for i,n in zip(Rzeta,range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[-2:,:],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[-2:,:],i) # -2: ???
        
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]

    zr=zr_[:,:,-1,:]

    # Add a level on top and bottom with no-gradient
    NO3,PO4=OGCM_Data['NO3'][:,:,-1,:],OGCM_Data['PO4'][:,:,-1,:]
    chl,phyt=OGCM_Data['chlorophyll'][:,:,-1,:],OGCM_Data['phytoplankton'][:,:,-1,:]
    oxygen=OGCM_Data['oxygen'][:,:,-1,:]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))

    PO4=np.hstack((np.expand_dims(PO4[:,0,:],axis=1)\
                    ,PO4,np.expand_dims(PO4[:,-1,:],axis=1)))

    chl=np.hstack((np.expand_dims(chl[:,0,:],axis=1)\
                    ,chl,np.expand_dims(chl[:,-1,:],axis=1)))

    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))

    oxygen=np.hstack((np.expand_dims(oxygen[:,0,:],axis=1)\
                    ,oxygen,np.expand_dims(oxygen[:,-1,:],axis=1)))



    print('!!! ztosigma_1d !!!')
    NO3_c,PO4_c=np.zeros_like(zr),np.zeros_like(zr)
    phyt_c,chl_c=np.zeros_like(zr),np.zeros_like(zr)
    oxygen_c=np.zeros_like(zr)

    for i,j,k,l,m,n in zip(NO3,PO4,chl,phyt,oxygen,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        PO4_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        chl_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));
        oxygen_c[n]=ru.ztosigma_1d(np.flip(m,axis=0),zr[n],np.flipud(Z));

    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_north'][:]=NO3_c
    ncI['PO4_north'][:]=PO4_c
    ncI['NH4_north'][:]=0
    ncI['TIC_north'][:]=2100
    ncI['oxygen_north'][:]=oxygen_c
    ncI['chlo_north'][:]=chl_c
    ncI['phyt_north'][:]=phyt_c
    ncI['zoop_north'][:]=phyt_c*0.3
    ncI['alkalinity_north'][:]=2350
   
    ncI['SdeC_north'][:]=0
    ncI['LdeC_north'][:]=0
    ncI['RdeC_north'][:]=0
    ncI['SdeN_north'][:]=0
    ncI['RdeN_north'][:]=0
    ncI['RdeN_north'][:]=0
    ncI.close()
    
    
if NSEW[1]:
    bry_lat_co=np.where( (latO[:]>=np.min(latG[:2,:])-latOmax) & (latO[:]<=np.max(latG[:2,:])+latOmax) )[0]
    bry_lon_co=lonO_co
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

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

        # tmp_data=np.squeeze(MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co])
        tmp_data=xr.open_mfdataset(OGCM_npth+'*.nc',parallel=Parallel)[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
        #tmp_data=MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co]  
        
        mask = np.isnan(tmp_data)
        DATA = np.ma.array(tmp_data, mask=mask)
        
        tmp_mask_=np.invert(DATA.mask)    
        data=np.zeros([len(Bry_time_num),len(depthO),2,lonG.shape[-1]])
    
        for t in tqdm(range(len(Bry_time_num))):
            for d in range(len(depthO)):
                # Interp mask
                tmp_mask=tmp_mask_[t,d]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                                DATA[t,d][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
                
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:2,:],latG[:2,:]) ,'cubic' )
                data[t,d]=data_re_.reshape(lonG[:2,:].shape) #.reshape(-1)
                
                tmp_var=data[t,d]
                
                if np.sum(~np.isnan(data[t,d]))==0:
                    data[t,d]=data[t,d-1]

                if np.sum(np.isnan(data[t,d]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t,d]=tmp_var
                    
        OGCM_Data[i]=data
        

    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    _,_,y_tmp,x_tmp=OGCM_Data['NO3'].shape

    Rzeta=np.zeros([len(TIMES_co),y_tmp,x_tmp])


    zr_= np.zeros([len(TIMES_co),MyVar['Layer_N'],Rzeta.shape[1],Rzeta.shape[-1]])
    zw = np.zeros([len(TIMES_co),MyVar['Layer_N']+1,Rzeta.shape[1],Rzeta.shape[-1]])

    
    for i,n in zip(Rzeta,range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[:2,:],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[:2,:],i) # -2: ???
        
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]

    zr=zr_[:,:,0,:]

    # Add a level on top and bottom with no-gradient
    NO3,PO4=OGCM_Data['NO3'][:,:,0,:],OGCM_Data['PO4'][:,:,0,:]
    chl,phyt=OGCM_Data['chlorophyll'][:,:,0,:],OGCM_Data['phytoplankton'][:,:,0,:]
    oxygen=OGCM_Data['oxygen'][:,:,0,:]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))

    PO4=np.hstack((np.expand_dims(PO4[:,0,:],axis=1)\
                    ,PO4,np.expand_dims(PO4[:,-1,:],axis=1)))

    chl=np.hstack((np.expand_dims(chl[:,0,:],axis=1)\
                    ,chl,np.expand_dims(chl[:,-1,:],axis=1)))

    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))

    oxygen=np.hstack((np.expand_dims(oxygen[:,0,:],axis=1)\
                    ,oxygen,np.expand_dims(oxygen[:,-1,:],axis=1)))

    print('!!! ztosigma_1d !!!')
    NO3_c,PO4_c=np.zeros_like(zr),np.zeros_like(zr)
    phyt_c,chl_c=np.zeros_like(zr),np.zeros_like(zr)
    oxygen_c=np.zeros_like(zr)

    for i,j,k,l,m,n in zip(NO3,PO4,chl,phyt,oxygen,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        PO4_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        chl_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));
        oxygen_c[n]=ru.ztosigma_1d(np.flip(m,axis=0),zr[n],np.flipud(Z));

    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_south'][:]=NO3_c
    ncI['PO4_south'][:]=PO4_c
    ncI['NH4_south'][:]=0
    ncI['TIC_south'][:]=2100
    ncI['oxygen_south'][:]=oxygen_c
    ncI['chlo_south'][:]=chl_c
    ncI['phyt_south'][:]=phyt_c
    ncI['zoop_south'][:]=phyt_c*0.3
    ncI['alkalinity_south'][:]=2350
   
    ncI['SdeC_south'][:]=0
    ncI['LdeC_south'][:]=0
    ncI['RdeC_south'][:]=0
    ncI['SdeN_south'][:]=0
    ncI['RdeN_south'][:]=0
    ncI['RdeN_south'][:]=0

    ncI.close()




if NSEW[2]:
    bry_lat_co=latO_co
    bry_lon_co=np.where( (lonO[:]>=np.min(lonG[:,-2:])-lonOmax) & (lonO[:]<=np.max(lonG[:,-2:])+lonOmax) )[0]
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

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

        # tmp_data=np.squeeze(MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co])
        tmp_data=xr.open_mfdataset(OGCM_npth+'*.nc',parallel=Parallel)[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
        #tmp_data=MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co]  
        
        mask = np.isnan(tmp_data)
        DATA = np.ma.array(tmp_data, mask=mask)
        
        tmp_mask_=np.invert(DATA.mask)    
        data=np.zeros([len(Bry_time_num),len(depthO),lonG.shape[0],2])
    
        for t in tqdm(range(len(Bry_time_num))):
            for d in range(len(depthO)):
                # Interp mask
                tmp_mask=tmp_mask_[t,d]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                                DATA[t,d][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
                
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:,-2:],latG[:,-2:]) ,'cubic' )
                data[t,d]=data_re_.reshape(lonG[:,-2:].shape) #.reshape(-1)
                
                tmp_var=data[t,d]
                
                if np.sum(~np.isnan(data[t,d]))==0:
                    data[t,d]=data[t,d-1]

                if np.sum(np.isnan(data[t,d]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t,d]=tmp_var
                    
        OGCM_Data[i]=data
        
    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    _,_,y_tmp,x_tmp=OGCM_Data['NO3'].shape

    Rzeta=np.zeros([len(TIMES_co),y_tmp,x_tmp])


    zr_= np.zeros([len(TIMES_co),MyVar['Layer_N'],Rzeta.shape[1],Rzeta.shape[-1]])
    zw = np.zeros([len(TIMES_co),MyVar['Layer_N']+1,Rzeta.shape[1],Rzeta.shape[-1]])

    
    for i,n in zip(Rzeta,range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[:,-2:],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[:,-2:],i) # -2: ???
        
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]

    zr=zr_[:,:,:,-1]

    # Add a level on top and bottom with no-gradient
    NO3,PO4=OGCM_Data['NO3'][:,:,:,-1],OGCM_Data['PO4'][:,:,:,-1]
    chl,phyt=OGCM_Data['chlorophyll'][:,:,:,-1],OGCM_Data['phytoplankton'][:,:,:,-1]
    oxygen=OGCM_Data['oxygen'][:,:,:,-1]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))

    PO4=np.hstack((np.expand_dims(PO4[:,0,:],axis=1)\
                    ,PO4,np.expand_dims(PO4[:,-1,:],axis=1)))

    chl=np.hstack((np.expand_dims(chl[:,0,:],axis=1)\
                    ,chl,np.expand_dims(chl[:,-1,:],axis=1)))

    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))

    oxygen=np.hstack((np.expand_dims(oxygen[:,0,:],axis=1)\
                    ,oxygen,np.expand_dims(oxygen[:,-1,:],axis=1)))

    print('!!! ztosigma_1d !!!')
    NO3_c,PO4_c=np.zeros_like(zr),np.zeros_like(zr)
    phyt_c,chl_c=np.zeros_like(zr),np.zeros_like(zr)
    oxygen_c=np.zeros_like(zr)

    for i,j,k,l,m,n in zip(NO3,PO4,chl,phyt,oxygen,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        PO4_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        chl_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));
        oxygen_c[n]=ru.ztosigma_1d(np.flip(m,axis=0),zr[n],np.flipud(Z));


    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_east'][:]=NO3_c
    ncI['PO4_east'][:]=PO4_c
    ncI['NH4_east'][:]=0
    ncI['TIC_east'][:]=2100
    ncI['oxygen_east'][:]=oxygen_c
    ncI['chlo_east'][:]=chl_c
    ncI['phyt_east'][:]=phyt_c
    ncI['zoop_east'][:]=phyt_c*0.3
    ncI['alkalinity_east'][:]=2350
   
    ncI['SdeC_east'][:]=0
    ncI['LdeC_east'][:]=0
    ncI['RdeC_east'][:]=0
    ncI['SdeN_east'][:]=0
    ncI['RdeN_east'][:]=0
    ncI['RdeN_east'][:]=0
    ncI.close()



    
    
 
if NSEW[3]:
    bry_lat_co=latO_co
    bry_lon_co=np.where( (lonO[:]>=np.min(lonG[:,:2])-lonOmax) & (lonO[:]<=np.max(lonG[:,:2])+lonOmax) )[0]
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

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

        # tmp_data=np.squeeze(MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co])
        tmp_data=xr.open_mfdataset(OGCM_npth+'*.nc',parallel=Parallel)[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
        #tmp_data=MFDataset(OGCM_npth+'*.nc')[OGCMVar[i]][TIMES_co,:,bry_lat_co,bry_lon_co]  
        
        mask = np.isnan(tmp_data)
        DATA = np.ma.array(tmp_data, mask=mask)
        
        tmp_mask_=np.invert(DATA.mask)    
        data=np.zeros([len(Bry_time_num),len(depthO),lonG.shape[0],2])
    
        for t in tqdm(range(len(Bry_time_num))):
            for d in range(len(depthO)):
                # Interp mask
                tmp_mask=tmp_mask_[t,d]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                                DATA[t,d][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
                
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:,:2],latG[:,:2]) ,'cubic' )
                data[t,d]=data_re_.reshape(lonG[:,:2].shape) #.reshape(-1)
                
                tmp_var=data[t,d]
                
                if np.sum(~np.isnan(data[t,d]))==0:
                    data[t,d]=data[t,d-1]

                if np.sum(np.isnan(data[t,d]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t,d]=tmp_var
                    
        OGCM_Data[i]=data
        
    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    _,_,y_tmp,x_tmp=OGCM_Data['NO3'].shape

    Rzeta=np.zeros([len(TIMES_co),y_tmp,x_tmp])


    zr_= np.zeros([len(TIMES_co),MyVar['Layer_N'],Rzeta.shape[1],Rzeta.shape[-1]])
    zw = np.zeros([len(TIMES_co),MyVar['Layer_N']+1,Rzeta.shape[1],Rzeta.shape[-1]])

    
    for i,n in zip(Rzeta,range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[:,:2],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[:,:2],i) # -2: ???
        
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]

    zr=zr_[:,:,:,0]

    # Add a level on top and bottom with no-gradient
    NO3,PO4=OGCM_Data['NO3'][:,:,:,0],OGCM_Data['PO4'][:,:,:,0]
    chl,phyt=OGCM_Data['chlorophyll'][:,:,:,0],OGCM_Data['phytoplankton'][:,:,:,0]
    oxygen=OGCM_Data['oxygen'][:,:,:,0]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))

    PO4=np.hstack((np.expand_dims(PO4[:,0,:],axis=1)\
                    ,PO4,np.expand_dims(PO4[:,-1,:],axis=1)))

    chl=np.hstack((np.expand_dims(chl[:,0,:],axis=1)\
                    ,chl,np.expand_dims(chl[:,-1,:],axis=1)))

    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))

    oxygen=np.hstack((np.expand_dims(oxygen[:,0,:],axis=1)\
                    ,oxygen,np.expand_dims(oxygen[:,-1,:],axis=1)))

    print('!!! ztosigma_1d !!!')
    NO3_c,PO4_c=np.zeros_like(zr),np.zeros_like(zr)
    phyt_c,chl_c=np.zeros_like(zr),np.zeros_like(zr)
    oxygen_c=np.zeros_like(zr)

    for i,j,k,l,m,n in zip(NO3,PO4,chl,phyt,oxygen,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        PO4_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        chl_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));
        oxygen_c[n]=ru.ztosigma_1d(np.flip(m,axis=0),zr[n],np.flipud(Z));


    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_west'][:]=NO3_c
    ncI['PO4_west'][:]=PO4_c
    ncI['NH4_west'][:]=0
    ncI['TIC_west'][:]=2100
    ncI['oxygen_west'][:]=oxygen_c
    ncI['chlo_west'][:]=chl_c
    ncI['phyt_west'][:]=phyt_c
    ncI['zoop_west'][:]=phyt_c*0.3
    ncI['alkalinity_west'][:]=2350
   
    ncI['SdeC_west'][:]=0
    ncI['LdeC_west'][:]=0
    ncI['RdeC_west'][:]=0
    ncI['SdeN_west'][:]=0
    ncI['RdeN_west'][:]=0
    ncI['RdeN_west'][:]=0
    ncI.close()

   
    
    
