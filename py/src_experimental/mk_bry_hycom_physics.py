
PKG_path = '/home/shjo/ROMS/romsforge/py/dev_individual/' # Location of JNUROMS directory
import sys 
sys.path.append(PKG_path)
import libs.ROMS_utils01 as ru
import libs.ROMS_utils02 as ru2
from libs.create_B import createB_NPZD, createB
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

My_Bry='/data/share/DATA/ROMS_INPUTS/bry/roms_bry_fennel_15km_smooth_v2.nc' # Initial file name (to create)
My_Grd='/data/share/DATA/ROMS_INPUTS/grd/roms_grd_fennel_15km_smooth_v2.nc' # Grd name

Parallel=False
#-- Define OGCM path ----------------------------------------------------------
ncdir='/data/share/DATA/RAW/00utc/'
sshNC=ncdir+'HYCOM_'
tempNC=ncdir+'HYCOM_'
saltNC=ncdir+'HYCOM_'
uNC=ncdir+'HYCOM_'
vNC=ncdir+'HYCOM_'

NSEW=[True,True,True,False] # N S E W
#NSEW=[False,False,False,True] # N S E W

Bry_title='ROMS-Fennel'

# OGCM Variables name
OGCMVar={'lon_rho':'lon','lat_rho':'lat','depth':'depth','time':'time',\
         'lon_u':'lon','lat_u':'lon','lon_v':'lon','lat_v':'lat',
         'temp':'water_temp','salt':'salinity','u':'water_u','v':'water_v','zeta':'surf_el'}

conserv=1

OGCMS=[ncdir+i for i in os.listdir(ncdir) if (i.endswith('.nc') & i.startswith('HYCOM_') ) ]

# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]

# My Variables
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

t_rng = ['2024-12-30 00:00', '2025-06-01 23:00']
My_time_ref = 'days since 2000-01-01 00:00:00'
TIME_UNIT = OGCM_TIMES.units

# Convert OGCM time to datetime (e.g. HYCOM time is in "seconds since 1970-01-01")
OGCM_times = num2date(OGCM_TIMES[:], TIME_UNIT)

# Define target datetime range
Tst = dt.datetime.strptime(t_rng[0], "%Y-%m-%d %H:%M")
Ted = dt.datetime.strptime(t_rng[1], "%Y-%m-%d %H:%M")

# Get the time indices in range
TIMES_co = np.where((OGCM_times >= Tst) & (OGCM_times <= Ted))[0]

# Convert to new reference time base
OGCM_days = date2num(OGCM_times[TIMES_co], TIME_UNIT)  # time in original unit
Bry_time_num = date2num(OGCM_times[TIMES_co], My_time_ref)  # re-referenced time
Bry_time_time = num2date(Bry_time_num, My_time_ref)

thO=ncO[OGCMVar['depth']].shape[0]
ncO.close()

print(TIMES_co)

atG,onG=lonG.shape
#cosa_,sina_=np.cos(angle)[-2:],np.sin(angle)[-2:] #NORTHERN BRY
#cosa=np.tile( np.tile(cosa_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
#sina=np.tile( np.tile(sina_,(thO,1,1)), (len(Bry_time_num),1,1,1) )

createB(My_Bry,topo,mask,MyVar,Bry_time_num,My_time_ref,NSEW,'Fennel','NETCDF3_64BIT_OFFSET')
#print("!!!=== #create ===!!!")

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
    print('--- Northern boundary ---')
    bry_lat_co=np.where( (latO[:]>=np.min(latG[-2:,:])-latOmax) & (latO[:]<=np.max(latG[-2:,:])+latOmax) )[0]
    bry_lon_co=lonO_co
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

    OGCM_Data={} #,OGCM_Mask={}

    for i in ['temp','u','v','salt','zeta','ubar','vbar']:
    #for i in ['temp']:
        print('!!! Data processing : '+i+' !!!')
        
        if (i=='zeta') or (i=='ubar') or (i=='vbar'):
            data=np.zeros([len(Bry_time_num),2,lonG.shape[-1]])

            if i=='zeta':
                # DATA=np.squeeze(MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co])
                DATA=xr.open_mfdataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co].squeeze().values
                #DATA=MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co]  

                mask = np.isnan(DATA)
                DATA = np.ma.array(DATA, mask=mask)
                tmp_mask_=np.invert(DATA.mask)
                
            elif i=='ubar':
                # tmp_u=np.squeeze(MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co])
                tmp_u=xr.open_mfdataset(uNC+'*.nc',parallel=Parallel)[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
                #tmp_u=MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co]  

                mask = np.isnan(tmp_u)
                tmp_u = np.ma.array(tmp_u, mask=mask)
                tmp_mask_=np.invert(tmp_u.mask)[:,:,:,:]
                
                tmp_u[tmp_u.mask]=0
                
                du=np.zeros([tmp_u.shape[0],tmp_u.shape[2],tmp_u.shape[3]])
                zu=np.zeros_like(du)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    du=du+dz[n]*tmp_u[:,n,:,:].data
                    zu=zu+dz[n]*tmp_mask_[:,n,:,:]
                DATA=du/zu
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            elif i=='vbar':
                #tmp_v=np.squeeze(MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co])
                #tmp_v=MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co]  
                tmp_v=xr.open_mfdataset(vNC+'*.nc',parallel=Parallel)[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values

                mask = np.isnan(tmp_v)
                tmp_v = np.ma.array(tmp_v, mask=mask)
                
                tmp_mask_=np.invert(tmp_v.mask)[:,:,:,:]
                
                tmp_v[tmp_v.mask]=0
                
                dv=np.zeros([tmp_v.shape[0],tmp_v.shape[2],tmp_v.shape[3]])
                zv=np.zeros_like(dv)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    dv=dv+dz[n]*tmp_v[:,n,:,:].data
                    zv=zv+dz[n]*tmp_mask_[:,n,:,:]
                DATA=dv/zv
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            for t in tqdm(range(len(Bry_time_num))):
                tmp_mask=tmp_mask_[t]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                                DATA[t][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
        
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[-2:,:].flatten(),latG[-2:,:].flatten()) ,'cubic' )
                data[t]=data_re_.reshape(lonG[-2:,:].shape)
                
                tmp_var=data[t]

                if np.sum(np.isnan(data[t]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t]=tmp_var
            OGCM_Data[i]=data
        else:    
            if i=='u' :
                OGCM_npth=uNC;
            elif i=='v':
                OGCM_npth=vNC;
            elif i=='temp':
                OGCM_npth=tempNC;
            elif i=='salt':
                OGCM_npth=saltNC;

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
    cosa_,sina_=np.cos(angle)[-2:,:],np.sin(angle)[-2:,:] #NORTHERN BRY
    cosa=np.tile( np.tile(cosa_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
    sina=np.tile( np.tile(sina_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
    #Process 2D vectors
    #modi

    ubar_north= ru2.rho2u_2d(OGCM_Data['ubar'][:,-1,:]*cosa[:,0,-1,:]+OGCM_Data['vbar'][:,-1,:]*sina[:,0,-1,:])
    vbar_north= OGCM_Data['vbar'][:,:,:]*cosa[:,0,:,:]-OGCM_Data['ubar'][:,:,:]*sina[:,0,:,:]
    vbar_north= (vbar_north[:,0,:]+vbar_north[:,-1,:])/2

    u=ru2.rho2u_3d(OGCM_Data['u'][:,:,-1,:]*cosa[:,:,-1,:]+OGCM_Data['v'][:,:,-1,:]*sina[:,:,-1,:])
    v=OGCM_Data['v'][:,:,:,:]*cosa[:,:,:,:]-OGCM_Data['u'][:,:,:,:]*sina[:,:,:,:]
    v=(v[:,:,0,:]+v[:,:,1,:])/2

    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    Rzeta=OGCM_Data['zeta'][:,-1,:] # -1 for northern BRY

    zr_=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N'],OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])
    zw=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N']+1,OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])

    for i,n in zip(OGCM_Data['zeta'],range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[-2:,:],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[-2:,:],i) # -2: ???

    zu= ru2.rho2u_3d(zr_[:,:,-1,:])
    zv=(zr_[:,:,0,:]+zr_[:,:,-1,:])/2

    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]

    dzu=ru2.rho2u_3d(dzr[:,:,-1,:])
    dzv=(dzr[:,:,0,:]+dzr[:,:,-1,:])/2
    zr=zr_[:,:,-1,:]

    # Add a level on top and bottom with no-gradient
    temp,salt=OGCM_Data['temp'][:,:,-1,:],OGCM_Data['salt'][:,:,-1,:]

    u1=np.hstack((np.expand_dims(u[:,0,:],axis=1)\
                    ,u,np.expand_dims(u[:,-1,:],axis=1)))
    v1=np.hstack((np.expand_dims(v[:,0,:],axis=1)\
                    ,v,np.expand_dims(v[:,-1,:],axis=1)))
    temp=np.hstack((np.expand_dims(temp[:,0,:],axis=1)\
                    ,temp,np.expand_dims(temp[:,-1,:],axis=1)))
    salt=np.hstack((np.expand_dims(salt[:,0,:],axis=1)\
                    ,salt,np.expand_dims(salt[:,-1,:],axis=1)))
    print('!!! ztosigma_1d !!!')
    u_c_,v_c_=np.zeros_like(zu),np.zeros_like(zv)
    temp_c,salt_c=np.zeros_like(zr),np.zeros_like(zr)
    for i,j,k,l,n in zip(u1,v1,temp,salt,range(zu.shape[0])): 
        u_c_[n]   =ru.ztosigma_1d(np.flip(i,axis=0),zu[n],np.flipud(Z));
        v_c_[n]   =ru.ztosigma_1d(np.flip(j,axis=0),zv[n],np.flipud(Z));
        temp_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        salt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));

    # Conservation
    if conserv==1:
        u_c= deepcopy(u_c_)  
        v_c= deepcopy(v_c_) 

        tmpu=np.sum(u_c_*dzu,axis=1)/np.sum(dzu,axis=1)
        tmpv=np.sum(v_c_*dzv,axis=1)/np.sum(dzv,axis=1)

        for i in range(dzu.shape[1]):
            u_c[:,i]=u_c[:,i,:] - tmpu +ubar_north
            v_c[:,i]=v_c[:,i,:] - tmpv +vbar_north

    # Barotropic velocities2
    ubar_north_=np.sum(u_c*dzu,axis=1)/np.sum(dzu,axis=1)
    vbar_north_=np.sum(v_c*dzv,axis=1)/np.sum(dzv,axis=1)

    ncI=Dataset(My_Bry,mode='a')
    ncI['zeta_north'][:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncI['temp_north'][:]=temp_c
    ncI['salt_north'][:]=salt_c
    ncI['u_north'][:]=u_c
    ncI['v_north'][:]=v_c
    ncI['ubar_north'][:]=ubar_north_
    ncI['vbar_north'][:]=vbar_north_
    ncI.close()

def new_func(data, t, tmp_var, d):
    data[t,d]=tmp_var

if NSEW[1]: # Southern bry
    print('--- Southern boundary ---')
    
    bry_lat_co=np.where( (latO[:]>=np.min(latG[:2,:])-latOmax) & (latO[:]<=np.max(latG[:2,:])+latOmax) )[0]
    bry_lon_co=lonO_co
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])

    OGCM_Data={}#,OGCM_Mask={}
    
    for i in ['temp','u','v','salt','zeta','ubar','vbar']:
    #for i in ['temp']:
        print('!!! Data processing : '+i+' !!!')
        
        if (i=='zeta') or (i=='ubar') or (i=='vbar'):
            data=np.zeros([len(Bry_time_num),2,lonG.shape[-1]])
            if i=='zeta':
                # DATA=np.squeeze(MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co])
                DATA=xr.open_mfdataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co].squeeze().values
                # DATA=MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co]  
    
                mask = np.isnan(DATA)
                DATA = np.ma.array(DATA, mask=mask)
                tmp_mask_=np.invert(DATA.mask)
                
            elif i=='ubar':
                # tmp_u=np.squeeze(MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co])
                tmp_u=xr.open_mfdataset(uNC+'*.nc',parallel=Parallel)[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
                # tmp_u=MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co]  
    
                mask = np.isnan(tmp_u)
                tmp_u = np.ma.array(tmp_u, mask=mask)
                tmp_mask_=np.invert(tmp_u.mask)[:,:,:,:]
                
                tmp_u[tmp_u.mask]=0
                
                du=np.zeros([tmp_u.shape[0],tmp_u.shape[2],tmp_u.shape[3]])
                zu=np.zeros_like(du)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    du=du+dz[n]*tmp_u[:,n,:,:].data
                    zu=zu+dz[n]*tmp_mask_[:,n,:,:]
                DATA=du/zu
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            elif i=='vbar':
                # tmp_v=np.squeeze(MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co])
                #tmp_v=MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co]  
                tmp_v=xr.open_mfdataset(vNC+'*.nc',parallel=Parallel)[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values

                mask = np.isnan(tmp_v)
                tmp_v = np.ma.array(tmp_v, mask=mask)
                
                tmp_mask_=np.invert(tmp_v.mask)[:,:,:,:]
                
                tmp_v[tmp_v.mask]=0
                
                dv=np.zeros([tmp_v.shape[0],tmp_v.shape[2],tmp_v.shape[3]])
                zv=np.zeros_like(dv)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    dv=dv+dz[n]*tmp_v[:,n,:,:].data
                    zv=zv+dz[n]*tmp_mask_[:,n,:,:]
                DATA=dv/zv
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            for t in tqdm(range(len(Bry_time_num))):
                tmp_mask=tmp_mask_[t]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                              DATA[t][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
        
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:2,:].flatten(),latG[:2,:].flatten()) ,'cubic' )
                data[t]=data_re_.reshape(lonG[:2,:].shape)
                
                tmp_var=data[t]
    
                if np.sum(np.isnan(data[t]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t]=tmp_var
            OGCM_Data[i]=data
        else:
            
            if i=='u' :
                OGCM_npth=uNC;
            elif i=='v':
                OGCM_npth=vNC;
            elif i=='temp':
                OGCM_npth=tempNC;
            elif i=='salt':
                OGCM_npth=saltNC;
                
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
                        new_func(data, t, tmp_var, d)
                        
            OGCM_Data[i]=data

    cosa_,sina_=np.cos(angle)[:2,:],np.sin(angle)[:2,:] #SOTHERN BRY
    cosa=np.tile( np.tile(cosa_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
    sina=np.tile( np.tile(sina_,(thO,1,1)), (len(Bry_time_num),1,1,1) )


    ubar_south= ru2.rho2u_2d(OGCM_Data['ubar'][:,0,:]*cosa[:,0,0,:]+OGCM_Data['vbar'][:,0,:]*sina[:,0,0,:])
    vbar_south= OGCM_Data['vbar'][:,:,:]*cosa[:,0,:,:]-OGCM_Data['ubar'][:,:,:]*sina[:,0,:,:]
    vbar_south= (vbar_south[:,1,:]+vbar_south[:,0,:])/2
    
    u=ru2.rho2u_3d(OGCM_Data['u'][:,:,0,:]*cosa[:,:,0,:]+OGCM_Data['v'][:,:,0,:]*sina[:,:,0,:])
    v=OGCM_Data['v'][:,:,:,:]*cosa[:,:,:,:]-OGCM_Data['u'][:,:,:,:]*sina[:,:,:,:]
    v=(v[:,:,0,:]+v[:,:,-1,:])/2
    
    
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    Rzeta=OGCM_Data['zeta'][:,0,:] # -1 for northern BRY
    zr_=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N'],OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])
    zw=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N']+1,OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])

    for i,n in zip(OGCM_Data['zeta'],range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[:2,:],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[:2,:],i) # -2: ???
    zu= ru2.rho2u_3d(zr_[:,:,0,:])
    zv=(zr_[:,:,0,:]+zr_[:,:,1,:])/2
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]
    dzu=ru2.rho2u_3d(dzr[:,:,0,:])
    dzv=(dzr[:,:,0,:]+dzr[:,:,1,:])/2
    zr=zr_[:,:,0,:]

    # Add a level on top and bottom with no-gradient
    temp,salt=OGCM_Data['temp'][:,:,0,:],OGCM_Data['salt'][:,:,0,:]

    u1=np.hstack((np.expand_dims(u[:,0,:],axis=1)\
                    ,u,np.expand_dims(u[:,-1,:],axis=1)))
    v1=np.hstack((np.expand_dims(v[:,0,:],axis=1)\
                    ,v,np.expand_dims(v[:,-1,:],axis=1)))
    temp=np.hstack((np.expand_dims(temp[:,0,:],axis=1)\
                    ,temp,np.expand_dims(temp[:,-1,:],axis=1)))
    salt=np.hstack((np.expand_dims(salt[:,0,:],axis=1)\
                    ,salt,np.expand_dims(salt[:,-1,:],axis=1)))
    print('!!! ztosigma_1d !!!')
    u_c_,v_c_=np.zeros_like(zu),np.zeros_like(zv)
    temp_c,salt_c=np.zeros_like(zr),np.zeros_like(zr)
    for i,j,k,l,n in zip(u1,v1,temp,salt,range(zu.shape[0])): 
        u_c_[n]   =ru.ztosigma_1d(np.flip(i,axis=0),zu[n],np.flipud(Z));
        v_c_[n]   =ru.ztosigma_1d(np.flip(j,axis=0),zv[n],np.flipud(Z));
        temp_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        salt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));

    # Conservation
    if conserv==1:
        u_c= deepcopy(u_c_)  
        v_c= deepcopy(v_c_) 

        tmpu=np.sum(u_c_*dzu,axis=1)/np.sum(dzu,axis=1)
        tmpv=np.sum(v_c_*dzv,axis=1)/np.sum(dzv,axis=1)

        for i in range(dzu.shape[1]):
            u_c[:,i]=u_c[:,i,:] - tmpu +ubar_south
            v_c[:,i]=v_c[:,i,:] - tmpv +vbar_south

    # Barotropic velocities2
    ubar_south_=np.sum(u_c*dzu,axis=1)/np.sum(dzu,axis=1)
    vbar_south_=np.sum(v_c*dzv,axis=1)/np.sum(dzv,axis=1)

    ncI=Dataset(My_Bry,mode='a')
    ncI['zeta_south'][:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncI['temp_south'][:]=temp_c
    ncI['salt_south'][:]=salt_c
    ncI['u_south'][:]=u_c
    ncI['v_south'][:]=v_c
    ncI['ubar_south'][:]=ubar_south_
    ncI['vbar_south'][:]=vbar_south_
    ncI.close()

if NSEW[2]: # Eastern bry
    print('--- Eastern boundary ---')
      
    bry_lat_co=latO_co
    bry_lon_co=np.where( (lonO[:]>=np.min(lonG[:,-2:])-lonOmax) & (lonO[:]<=np.max(lonG[:,-2:])+lonOmax) )[0]
    
    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])
    
    OGCM_Data={}#,OGCM_Mask={}
    
    for i in ['temp','u','v','salt','zeta','ubar','vbar']:
    #for i in ['temp']:
        print('!!! Data processing : '+i+' !!!')
        
        if (i=='zeta') or (i=='ubar') or (i=='vbar'):
            data=np.zeros([len(Bry_time_num),lonG.shape[0],2])
            if i=='zeta':
                # DATA=np.squeeze(MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co])
                DATA=xr.open_mfdataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co].squeeze().values
                #DATA=MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co]  
    
                mask = np.isnan(DATA)
                DATA = np.ma.array(DATA, mask=mask)
                tmp_mask_=np.invert(DATA.mask)
                
            elif i=='ubar':
                # tmp_u=np.squeeze(MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co])
                tmp_u=xr.open_mfdataset(uNC+'*.nc',parallel=Parallel)[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
                #tmp_u=MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co]  
    
                mask = np.isnan(tmp_u)
                tmp_u = np.ma.array(tmp_u, mask=mask)
                tmp_mask_=np.invert(tmp_u.mask)[:,:,:,:]
                
                tmp_u[tmp_u.mask]=0
                
                du=np.zeros([tmp_u.shape[0],tmp_u.shape[2],tmp_u.shape[3]])
                zu=np.zeros_like(du)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    du=du+dz[n]*tmp_u[:,n,:,:].data
                    zu=zu+dz[n]*tmp_mask_[:,n,:,:]
                DATA=du/zu
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            elif i=='vbar':
                # tmp_v=np.squeeze(MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co])
                #tmp_v=MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co]  
                tmp_v=xr.open_mfdataset(vNC+'*.nc',parallel=Parallel)[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values

                mask = np.isnan(tmp_v)
                tmp_v = np.ma.array(tmp_v, mask=mask)
                
                tmp_mask_=np.invert(tmp_v.mask)[:,:,:,:]
                
                tmp_v[tmp_v.mask]=0
                
                dv=np.zeros([tmp_v.shape[0],tmp_v.shape[2],tmp_v.shape[3]])
                zv=np.zeros_like(dv)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    dv=dv+dz[n]*tmp_v[:,n,:,:].data
                    zv=zv+dz[n]*tmp_mask_[:,n,:,:]
                DATA=dv/zv
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            for t in tqdm(range(len(Bry_time_num))):
                tmp_mask=tmp_mask_[t]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                              DATA[t][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
        
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:,-2:].flatten(),latG[:,-2:].flatten()) ,'cubic' )
                data[t]=data_re_.reshape(lonG[:,-2:].shape)
                
                tmp_var=data[t]
    
                if np.sum(np.isnan(data[t]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t]=tmp_var
            OGCM_Data[i]=data
        else:
            
            if i=='u' :
                OGCM_npth=uNC;
            elif i=='v':
                OGCM_npth=vNC;
            elif i=='temp':
                OGCM_npth=tempNC;
            elif i=='salt':
                OGCM_npth=saltNC;
                
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

    cosa_,sina_=np.cos(angle)[:,-2:],np.sin(angle)[:,-2:] #EASTERN BRY
    cosa=np.tile( np.tile(cosa_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
    sina=np.tile( np.tile(sina_,(thO,1,1)), (len(Bry_time_num),1,1,1) )

    
    ubar_east= OGCM_Data['ubar'][:,:,:]*cosa[:,0,:,:]+OGCM_Data['vbar'][:,:,:]*sina[:,0,:,:]
    ubar_east= (ubar_east[:,:,0]+ubar_east[:,:,-1])/2
    vbar_east=ru2.rho2u_2d(OGCM_Data['vbar'][:,:,-1]*cosa[:,0,:,-1]-OGCM_Data['ubar'][:,:,-1]*sina[:,0,:,-1])
    
    u=OGCM_Data['u'][:,:,:,:]*cosa[:,:,:,:]+OGCM_Data['v'][:,:,:,:]*sina[:,:,:,:]
    u=(u[:,:,:,0]+u[:,:,:,-1])/2
    v=ru2.rho2u_3d(OGCM_Data['v'][:,:,:,-1]*cosa[:,:,:,-1]-OGCM_Data['u'][:,:,:,-1]*sina[:,:,:,-1])

    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    Rzeta=OGCM_Data['zeta'][:,:,-1] # -1 for eastern BRY

    zr_=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N'],OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])
    zw=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N']+1,OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])

    for i,n in zip(OGCM_Data['zeta'],range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[:,-2:],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[:,-2:],i) # -2: ???
    zu= (zr_[:,:,:,-1]+zr_[:,:,:,0])/2
    zv=ru2.rho2u_3d(zr_[:,:,:,-1])
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]
    dzu=(dzr[:,:,:,-1]+dzr[:,:,:,0])/2
    dzv=ru2.rho2u_3d(dzr[:,:,:,-1])
    zr=zr_[:,:,:,-1]

    # Add a level on top and bottom with no-gradient
    temp,salt=OGCM_Data['temp'][:,:,:,-1],OGCM_Data['salt'][:,:,:,-1]

    u1=np.hstack((np.expand_dims(u[:,0,:],axis=1)\
                    ,u,np.expand_dims(u[:,-1,:],axis=1)))
    v1=np.hstack((np.expand_dims(v[:,0,:],axis=1)\
                    ,v,np.expand_dims(v[:,-1,:],axis=1)))
    temp=np.hstack((np.expand_dims(temp[:,0,:],axis=1)\
                    ,temp,np.expand_dims(temp[:,-1,:],axis=1)))
    salt=np.hstack((np.expand_dims(salt[:,0,:],axis=1)\
                    ,salt,np.expand_dims(salt[:,-1,:],axis=1)))

    print('!!! ztosigma_1d !!!')
    u_c_,v_c_=np.zeros_like(zu),np.zeros_like(zv)
    temp_c,salt_c=np.zeros_like(zr),np.zeros_like(zr)
    for i,j,k,l,n in zip(u1,v1,temp,salt,range(zu.shape[0])): 
        u_c_[n]   =ru.ztosigma_1d(np.flip(i,axis=0),zu[n],np.flipud(Z));
        v_c_[n]   =ru.ztosigma_1d(np.flip(j,axis=0),zv[n],np.flipud(Z));
        temp_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        salt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));
    # =============================================================================
    # Conservation
    if conserv==1:
        u_c= deepcopy(u_c_)  
        v_c= deepcopy(v_c_) 

        tmpu=np.sum(u_c_*dzu,axis=1)/np.sum(dzu,axis=1)
        tmpv=np.sum(v_c_*dzv,axis=1)/np.sum(dzv,axis=1)

        for i in range(dzu.shape[1]):
            u_c[:,i]=u_c[:,i,:] - tmpu +ubar_east
            v_c[:,i]=v_c[:,i,:] - tmpv +vbar_east

    # Barotropic velocities2
    ubar_east_=np.sum(u_c*dzu,axis=1)/np.sum(dzu,axis=1)
    vbar_east_=np.sum(v_c*dzv,axis=1)/np.sum(dzv,axis=1)

    ncI=Dataset(My_Bry,mode='a')
    ncI['zeta_east'][:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncI['temp_east'][:]=temp_c
    ncI['salt_east'][:]=salt_c
    ncI['u_east'][:]=u_c
    ncI['v_east'][:]=v_c
    ncI['ubar_east'][:]=ubar_east_
    ncI['vbar_east'][:]=vbar_east_
    ncI.close()

if NSEW[3]: # Western bry
    print('--- Western boundary ---')

    # West (gridbuilder version)
    bry_lat_co=latO_co
    bry_lon_co=np.where( (lonO[:]>=np.min(lonG[:,:2])-lonOmax) & (lonO[:]<=np.max(lonG[:,:2])+lonOmax) )[0]

    lonO_s_m,latO_s_m=np.meshgrid(lonO[bry_lon_co],latO[bry_lat_co])


    OGCM_Data={}#,OGCM_Mask={}

    for i in ['temp','u','v','salt','zeta','ubar','vbar']:
    #for i in ['temp']:
        print('!!! Data processing : '+i+' !!!')
        
        if (i=='zeta') or (i=='ubar') or (i=='vbar'):
            data=np.zeros([len(Bry_time_num),lonG.shape[0],2])
            if i=='zeta':
                # DATA=np.squeeze(MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co])
                DATA=xr.open_mfdataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co].squeeze().values
                # DATA=MFDataset(sshNC+'*.nc')[OGCMVar[i]][TIMES_co,bry_lat_co,bry_lon_co]  

                mask = np.isnan(DATA)
                DATA = np.ma.array(DATA, mask=mask)
                tmp_mask_=np.invert(DATA.mask)
                
            elif i=='ubar':
                # tmp_u=np.squeeze(MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co])
                tmp_u=xr.open_mfdataset(uNC+'*.nc',parallel=Parallel)[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values
                # tmp_u=MFDataset(uNC+'*.nc')[OGCMVar['u']][TIMES_co,:,bry_lat_co,bry_lon_co]  

                mask = np.isnan(tmp_u)
                tmp_u = np.ma.array(tmp_u, mask=mask)
                tmp_mask_=np.invert(tmp_u.mask)[:,:,:,:]
                
                tmp_u[tmp_u.mask]=0
                
                du=np.zeros([tmp_u.shape[0],tmp_u.shape[2],tmp_u.shape[3]])
                zu=np.zeros_like(du)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    du=du+dz[n]*tmp_u[:,n,:,:].data
                    zu=zu+dz[n]*tmp_mask_[:,n,:,:]
                DATA=du/zu
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            elif i=='vbar':
                # tmp_v=np.squeeze(MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co])
                #tmp_v=MFDataset(vNC+'*.nc')[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co]  
                tmp_v=xr.open_mfdataset(vNC+'*.nc',parallel=Parallel)[OGCMVar['v']][TIMES_co,:,bry_lat_co,bry_lon_co].squeeze().values

                mask = np.isnan(tmp_v)
                tmp_v = np.ma.array(tmp_v, mask=mask)
                
                tmp_mask_=np.invert(tmp_v.mask)[:,:,:,:]
                
                tmp_v[tmp_v.mask]=0
                
                dv=np.zeros([tmp_v.shape[0],tmp_v.shape[2],tmp_v.shape[3]])
                zv=np.zeros_like(dv)
                dz=np.gradient(-depthO)
                for n in range(len(depthO)):
                    dv=dv+dz[n]*tmp_v[:,n,:,:].data
                    zv=zv+dz[n]*tmp_mask_[:,n,:,:]
                DATA=dv/zv
                # DATA[DATA==0]=np.nan
                tmp_mask_=tmp_mask_[:,0,:,:]
                
            for t in tqdm(range(len(Bry_time_num))):
                tmp_mask=tmp_mask_[t]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                              DATA[t][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
        
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:,:2].flatten(),latG[:,:2].flatten()) ,'cubic' )
                data[t]=data_re_.reshape(lonG[:,:2].shape)
                
                tmp_var=data[t]

                if np.sum(np.isnan(data[t]))!=0:
                    tmp_var[np.isnan(tmp_var)]=np.nanmean(tmp_var)
                    data[t]=tmp_var
            OGCM_Data[i]=data
        else:
            
            if i=='u' :
                OGCM_npth=uNC;
            elif i=='v':
                OGCM_npth=vNC;
            elif i=='temp':
                OGCM_npth=tempNC;
            elif i=='salt':
                OGCM_npth=saltNC;
                
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
    cosa_,sina_=np.cos(angle)[:,:2],np.sin(angle)[:,:2] #WESTERN BRY
    cosa=np.tile( np.tile(cosa_,(thO,1,1)), (len(Bry_time_num),1,1,1) )
    sina=np.tile( np.tile(sina_,(thO,1,1)), (len(Bry_time_num),1,1,1) )

    ubar_west= OGCM_Data['ubar'][:,:,:]*cosa[:,0,:,:]+OGCM_Data['vbar'][:,:,:]*sina[:,0,:,:]
    ubar_west=(ubar_west[:,:,0]+ubar_west[:,:,-1])/2
    vbar_west=ru2.rho2u_2d(OGCM_Data['vbar'][:,:,0]*cosa[:,0,:,0]-OGCM_Data['ubar'][:,:,0]*sina[:,0,:,0])
    u=OGCM_Data['u'][:,:,:,:]*cosa[:,:,:,:]+OGCM_Data['v'][:,:,:,:]*sina[:,:,:,:]
    u=(u[:,:,:,0]+u[:,:,:,-1])/2
    v=ru2.rho2u_3d(OGCM_Data['v'][:,:,:,0]*cosa[:,:,:,0]-OGCM_Data['u'][:,:,:,0]*sina[:,:,:,0])

    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

    Rzeta=OGCM_Data['zeta'][:,:,0] # -1 for Western BRY


    zr_=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N'],OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])
    zw=np.zeros([OGCM_Data['zeta'].shape[0],MyVar['Layer_N']+1,OGCM_Data['zeta'].shape[1],OGCM_Data['zeta'].shape[-1]])

    for i,n in zip(OGCM_Data['zeta'],range(Rzeta.shape[0])):
        zr_[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        1,topo[:,:2],i);  # -2: ???    
        zw[n,:,:,:]=ru.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
                    MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                        5,topo[:,:2],i) # -2: ???
    zu=(zr_[:,:,:,0]+zr_[:,:,:,-1])/2
    zv=ru2.rho2u_3d(zr_[:,:,:,0])
    dzr=zw[:,1:,:,:]-zw[:,:-1,:,:] # [t,depth,lat,lon]
    dzu=(dzr[:,:,:,0]+dzr[:,:,:,-1])/2
    dzv=ru2.rho2u_3d(dzr[:,:,:,0])
    zr=zr_[:,:,:,0]
    # Add a level on top and bottom with no-gradient
    temp,salt=OGCM_Data['temp'][:,:,:,0],OGCM_Data['salt'][:,:,:,0]

    u1=np.hstack((np.expand_dims(u[:,0,:],axis=1)\
                    ,u,np.expand_dims(u[:,-1,:],axis=1)))
    v1=np.hstack((np.expand_dims(v[:,0,:],axis=1)\
                    ,v,np.expand_dims(v[:,-1,:],axis=1)))
    temp=np.hstack((np.expand_dims(temp[:,0,:],axis=1)\
                    ,temp,np.expand_dims(temp[:,-1,:],axis=1)))
    salt=np.hstack((np.expand_dims(salt[:,0,:],axis=1)\
                    ,salt,np.expand_dims(salt[:,-1,:],axis=1)))

    print('!!! ztosigma_1d !!!')
    u_c_,v_c_=np.zeros_like(zu),np.zeros_like(zv)
    temp_c,salt_c=np.zeros_like(zr),np.zeros_like(zr)
    for i,j,k,l,n in zip(u1,v1,temp,salt,range(zu.shape[0])): 
        u_c_[n]   =ru.ztosigma_1d(np.flip(i,axis=0),zu[n],np.flipud(Z));
        v_c_[n]   =ru.ztosigma_1d(np.flip(j,axis=0),zv[n],np.flipud(Z));
        temp_c[n]=ru.ztosigma_1d(np.flip(k,axis=0),zr[n],np.flipud(Z));
        salt_c[n]=ru.ztosigma_1d(np.flip(l,axis=0),zr[n],np.flipud(Z));

    # Conservation
    if conserv==1:
        u_c= deepcopy(u_c_)  
        v_c= deepcopy(v_c_) 

        tmpu=np.sum(u_c_*dzu,axis=1)/np.sum(dzu,axis=1)
        tmpv=np.sum(v_c_*dzv,axis=1)/np.sum(dzv,axis=1)

        for i in range(dzu.shape[1]):
            u_c[:,i]=u_c[:,i,:] - tmpu +ubar_west
            v_c[:,i]=v_c[:,i,:] - tmpv +vbar_west

    # Barotropic velocities2
    ubar_west_=np.sum(u_c*dzu,axis=1)/np.sum(dzu,axis=1)
    vbar_west_=np.sum(v_c*dzv,axis=1)/np.sum(dzv,axis=1)

    ncI=Dataset(My_Bry,mode='a')
    ncI['zeta_west'][:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncI['temp_west'][:]=temp_c
    ncI['salt_west'][:]=salt_c
    ncI['u_west'][:]=u_c
    ncI['v_west'][:]=v_c
    ncI['ubar_west'][:]=ubar_west_
    ncI['vbar_west'][:]=vbar_west_
    ncI.close()
