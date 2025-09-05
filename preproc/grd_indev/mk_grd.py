# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 12:48:19 2022

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""
import sys
sys.path.append('D:/Working_hub/OneDrive/base142/Factory/MantaROMS/src_d')
import JNU_CLASS2 as jr
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np

topo_path = 'D:/Working_hub/OneDrive/Sources/Data/Etopo/'
topo_name = topo_path+'etopo2.nc'
dl=1
lon_rng,lat_rng = [112,145], [10,30] 

# =============================================================================
# 
# =============================================================================

# Make isotropic grid
lonr=np.arange(lon_rng[0],lon_rng[-1],dl)

latr=[lat_rng[0]]
i=0;
while latr[i]<=lat_rng[-1]:
    i+=1
    latr.append(latr[i-1]+dl*np.cos(latr[i-1]*np.pi/180))

Lonr,Latr=np.meshgrid(lonr,latr)
Lonu,Lonv,Lonp=jr.rho2uvp(Lonr)
Latu,Latv,Latp=jr.rho2uvp(Latr)

M,L=Latp.shape
# =============================================================================
import JNU_CLASS2 as jr
# =============================================================================

# Compute the metrics
pm,pn,dndx,dmde=jr.get_metrics(Latu,Lonu,Latv,Lonv)
xr,yr=0*pm,0*pm
for i in range(L):
    xr[:,i+1]=xr[:,i]+2/(pm[:,i+1]+pm[:,i])
for j in range(M):
    yr[j+1,:]=yr[j,:]+2/(pn[j+1,:]+pn[j,:])
    
xu,xv,xp=jr.rho2uvp(xr)
yu,yv,yp=jr.rho2uvp(yr)
dx,dy=1/pm,1/pn
dxmax,dxmin=np.max(dx/1000),np.min(dx/1000)
dymax,dymin=np.max(dy/1000),np.min(dy/1000)
print(' Min dx='+str(dxmin)+' km - Max dx='+str(dxmax)+' km')
print(' Min dy='+str(dymin)+' km - Max dy='+str(dymax)+' km')

# =============================================================================
import JNU_CLASS2 as jr
# =============================================================================

# Angle between XI-axis and the direction to the East at RHO-points [radians]
angle=jr.get_angle(Latu,Lonu)
# Coriolis parameter
f=4*np.pi*np.sin(np.pi*Latr/180)/(24*3600);

# =============================================================================
# 
# =============================================================================

h=jr.add_topo(topo_name,lat_rng,lon_rng,pm,pn,Lonr,Latr)
maskr=h>-30
# Add topo
# maskr=process_mask(maskr);
masku,maskv,maskp=jr.uvp_mask(maskr);



import matplotlib.pyplot as plt
plt.pcolor(h)
plt.pcolor(maskr)






