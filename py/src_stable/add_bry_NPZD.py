
PKG_path = 'C:/Users/ust21/shjo/projects/myROMS/prc_tools/' # Location of JNUROMS directory
import sys 
sys.path.append(PKG_path)
import utils.ROMS_utils01 as ru
import utils.ROMS_utils02 as ru2
from utils.ncCreate import create_bry_ust,create_bry_NPZD
import xarray as xr
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import netCDF4 as nc4
import os
from scipy.interpolate import griddata
from copy import deepcopy
import datetime as dt
from tqdm import tqdm
import matplotlib.pyplot as plt

My_Bry='D:/shjo/ROMS_inputs/NWP15km/N36/ROMS_bry_15km_241231_250501_N36.nc' # Initial file name (to create)
My_Grd='D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc' # Grd name
 
Parallel=False
#-- Define OGCM path ----------------------------------------------------------
ncdir='D:/shjo/GLORYS/'
NO3NC=ncdir+'nut/mercatorbiomer4v2r1_global_mean_nut_'
phytNC=ncdir+'chl/mercatorbiomer4v2r1_global_mean_pft_'
zoopNC=ncdir+'zoo/mercatorbiomer4v2r1_global_mean_plankton_'
detrNC=ncdir+'HYCOM_'

NSEW=[True,True,True,False] # N S E W

Bry_title='test'

# OGCM Variables name
OGCMVar={'lon_rho':'longitude','lat_rho':'latitude','depth':'depth','time':'time',\
         'NO3':'no3','phyt':'phyc','zoop':'zooc','detr':'???'}

conserv=1

OGCMS=[ncdir+'nut/'+i for i in os.listdir(ncdir+'nut/') if (i.endswith('.nc') & i.startswith('mercatorbiomer4v2r1_global_mean_nut_') ) ]

# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]

MyVar={'Layer_N':36,'Vtransform':2,\
       'Vstretching':4,'Theta_s':6.5,\
           'Theta_b':1,'Tcline':400,'hmin':10}
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

t_rng = ['2024-12-31 00:00', '2025-05-01 23:00']
# t_rng = ['2024-12-31 00:00', '2025-04-01 23:00']
My_time_ref = 'days since 2000-01-01 00:00:00'
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
#cosa_,sina_=np.cos(angle)[-2:],np.sin(angle)[-2:] #NORTHERN BRY
#cosa=np.tile( np.tile(cosa_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
#sina=np.tile( np.tile(sina_,(thO,1,1)), (len(Bry_time_num),1,1,1) )

#create_bry_NPZD(My_Bry,mask,topo,MyVar,NSEW,OGCM_TIMES[:],TIME_UNIT,Bry_title,ncFormat='NETCDF4')


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

    for i in ['NO3','phyt','zoop']:
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
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[-2:,:],latG[-2:,:]) ,'cubic' )
                data[t,d]=data_re_.reshape(lonG[-2:,:].shape) #.reshape(-1)
                
                tmp_var=data[t,d]
                
                if np.sum(~np.isnan(data[t,d]))==0:
                    data[t,d]=data[t,d-1]

                if np.sum(np.isnan(data[t,d]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t,d]=tmp_var
                    
        OGCM_Data[i]=data
        
    OGCM_Data['detr']=OGCM_Data['phyt']*0.1

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
    NO3,phyt=OGCM_Data['NO3'][:,:,-1,:],OGCM_Data['phyt'][:,:,-1,:]
    zoop,detr=OGCM_Data['zoop'][:,:,-1,:],OGCM_Data['detr'][:,:,-1,:]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))
    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))
    zoop=np.hstack((np.expand_dims(zoop[:,0,:],axis=1)\
                    ,zoop,np.expand_dims(zoop[:,-1,:],axis=1)))
    detr=np.hstack((np.expand_dims(detr[:,0,:],axis=1)\
                    ,detr,np.expand_dims(detr[:,-1,:],axis=1)))
               
    print('!!! ztosigma_1d !!!')
    NO3_c,phyt_c=np.zeros_like(zr),np.zeros_like(zr)
    zoop_c,detr_c=np.zeros_like(zr),np.zeros_like(zr)

    for i,j,k,l,n in zip(NO3,phyt,zoop,detr,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        zoop_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        detr_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));

    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_north'][:]=NO3_c
    ncI['phyt_north'][:]=phyt_c
    ncI['zoop_north'][:]=zoop_c
    ncI['detritus_north'][:]=0.04
    ncI.close()
    
    
if NSEW[1]:
    bry_lat_co=np.where( (latO[:]>=np.min(latG[:2,:])-latOmax) & (latO[:]<=np.max(latG[:2,:])+latOmax) )[0]
    bry_lon_co=lonO_co
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

    OGCM_Data={}#,OGCM_Mask={}

    for i in ['NO3','phyt','zoop']:
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
        
    OGCM_Data['detr']=OGCM_Data['phyt']*0.1

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
    NO3,phyt=OGCM_Data['NO3'][:,:,0,:],OGCM_Data['phyt'][:,:,0,:]
    zoop,detr=OGCM_Data['zoop'][:,:,0,:],OGCM_Data['detr'][:,:,0,:]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))
    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))
    zoop=np.hstack((np.expand_dims(zoop[:,0,:],axis=1)\
                    ,zoop,np.expand_dims(zoop[:,-1,:],axis=1)))
    detr=np.hstack((np.expand_dims(detr[:,0,:],axis=1)\
                    ,detr,np.expand_dims(detr[:,-1,:],axis=1)))
               
    print('!!! ztosigma_1d !!!')
    NO3_c,phyt_c=np.zeros_like(zr),np.zeros_like(zr)
    zoop_c,detr_c=np.zeros_like(zr),np.zeros_like(zr)

    for i,j,k,l,n in zip(NO3,phyt,zoop,detr,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        zoop_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        detr_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));

    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_south'][:]=NO3_c
    ncI['phyt_south'][:]=phyt_c
    ncI['zoop_south'][:]=zoop_c
    ncI['detritus_south'][:]=0.04
    ncI.close()


if NSEW[2]:
    bry_lat_co=latO_co
    bry_lon_co=np.where( (lonO[:]>=np.min(lonG[:,-2:])-lonOmax) & (lonO[:]<=np.max(lonG[:,-2:])+lonOmax) )[0]
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

    OGCM_Data={}#,OGCM_Mask={}

    for i in ['NO3','phyt','zoop']:
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
        
    OGCM_Data['detr']=OGCM_Data['phyt']*0.1

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
    NO3,phyt=OGCM_Data['NO3'][:,:,:,-1],OGCM_Data['phyt'][:,:,:,-1]
    zoop,detr=OGCM_Data['zoop'][:,:,:,-1],OGCM_Data['detr'][:,:,:,-1]

    NO3=np.hstack((np.expand_dims(NO3[:,0,:],axis=1)\
                    ,NO3,np.expand_dims(NO3[:,-1,:],axis=1)))
    phyt=np.hstack((np.expand_dims(phyt[:,0,:],axis=1)\
                    ,phyt,np.expand_dims(phyt[:,-1,:],axis=1)))
    zoop=np.hstack((np.expand_dims(zoop[:,0,:],axis=1)\
                    ,zoop,np.expand_dims(zoop[:,-1,:],axis=1)))
    detr=np.hstack((np.expand_dims(detr[:,0,:],axis=1)\
                    ,detr,np.expand_dims(detr[:,-1,:],axis=1)))
               
    print('!!! ztosigma_1d !!!')
    NO3_c,phyt_c=np.zeros_like(zr),np.zeros_like(zr)
    zoop_c,detr_c=np.zeros_like(zr),np.zeros_like(zr)

    for i,j,k,l,n in zip(NO3,phyt,zoop,detr,range(zr.shape[0])): 
        NO3_c[n]=ru.ztosigma_1d(np.flip(i,axis=0),zr[n],np.flipud(Z));
        phyt_c[n]=ru.ztosigma_1d(np.flip(j,axis=0),zr[n],np.flipud(Z));
        zoop_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        detr_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));

    ncI=Dataset(My_Bry,mode='a')
    ncI['NO3_east'][:]=NO3_c
    ncI['phyt_east'][:]=phyt_c
    ncI['zoop_east'][:]=zoop_c
    ncI['detritus_east'][:]=0.04
    ncI.close()



    
    
    
    
    
