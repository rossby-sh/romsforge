# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 10:56:45 2025

@author: ust21
"""

import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import numpy as np
import sys
sys.path.append('C:/Users/ust21/shjo/others/will_delet/myROMS/prc_src/utils/')
from netCDF4 import Dataset, num2date
import ROMS_utils02 as ru
import cmocean

# Gpth='D:/shjo/ROMS_inputs/NWP12_grd_NWP4.nc'
# Gpth='D:/shjo/ROMS_inputs/roms_grd_fennel_5km_smooth_v2.nc'

# Gpth='D:/shjo/ROMS_inputs/NWP12_grd_NWP4.nc'

# Spth='D:/shjo/MCC/nifs_jun_srf3/OSTIA_SST_250601-250630.nc'
Spth='D:/shjo/ROMS_forcings/ROMS_F_ERA5_3hourly_250531-250701.nc'

import xarray as xr

ncX=xr.open_dataset(Spth).sst.loc[dict(sst_time=slice('2025-06-01','2025-06-30'))]
SST=ncX.mean(dim='sst_time')

title=''

wpth='D:/shjo/ROMS_forcings/'

# flist=[Spth+i for i in os.listdir(Spth) if i.endswith('.nc')]
# flist=[Spth+i for i in os.listdir(Spth) if i.startswith('NWP12_avg_')]
# flist=[Spth+i for i in os.listdir(Spth) if i.startswith('NIFS')]

# 0226 8457


def draw_roms_pcolor(lon2d, lat2d, var2d, timestamp=None,
                     varname='Variable', units='', log_scale=False,
                     clim=None, output_path=None, cmap='jet'):
    """
    ROMS 2D 결과를 pcolormesh로 시각화합니다.

    Parameters
    ----------
    lon2d, lat2d : 2D np.ndarray
        ROMS grid (rho-points 기준)
    var2d : 2D np.ndarray
        시각화할 변수 (예: chl[t, :, :])
    timestamp : datetime.datetime or str, optional
        제목에 사용할 날짜
    varname : str
        변수 이름
    units : str
        단위
    log_scale : bool
        log10 변환 여부
    clim : tuple of (vmin, vmax), optional
        컬러 스케일 범위
    output_path : str
        저장 경로. None이면 화면에 표시만 함.
    cmap : str
        사용할 컬러맵 (기본: 'jet')
    """
    # 값 전처리
    if log_scale:
        val_plot = np.where(var2d > 1e-20, np.log10(var2d), np.nan) if log_scale else var2d.copy()
    else:
        val_plot = np.where(var2d , var2d, np.nan) if log_scale else var2d.copy()

    # val_plot = np.where(var2d , np.log10(var2d), np.nan) if log_scale else var2d.copy()

    # 색상 범위
    if clim is None:
        vmin, vmax = -2, 1.3
    else:
        vmin, vmax = clim

    # val_plot[val_plot>=vmax-0.01]=vmax-0.01
    # val_plot[val_plot<=vmin+0.01]=vmin+0.01
    # print(np.nanmin(val_plot),np.nanmax(val_plot))

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    plt.rcParams.update({
        'font.size': 16,
        'axes.titlesize': 20,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 15,
        'figure.titlesize': 20
    })
    # 지도 설정
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([115, 150, 20, 49], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.')
    gl.top_labels = False
    gl.right_labels = False

    ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=[4,8,12,16,20,24,28,30])
    ax.clabel(ct, colors='k', fontsize=12)
    # pcolormesh
    pcm = ax.pcolormesh(
        lon2d, lat2d, val_plot,
        cmap=cmap, norm=norm, shading='auto',
        transform=ccrs.PlateCarree()
    )

    # 제목
    title = f'{varname}'
    if timestamp is not None:
        if hasattr(timestamp, 'strftime'):
            title += f' - {timestamp.strftime("%Y-%m-%d")}'
        else:
            title += f' - {timestamp}'
    ax.set_title(title, fontsize=16, pad=10)

    # colorbar
    cbar = plt.colorbar(pcm, ax=ax, orientation='vertical', pad=0.02, shrink=0.8)
    label = f'{"log(" if log_scale else ""}{varname}{")" if log_scale else ""} [{units}]'
    cbar.set_label(label, fontsize=14)

    # log tick (optional)
    if log_scale:
        tick_vals = [-2, -1.52, -1, -0.52, 0, 0.48, 1]
        tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10']
        cbar.set_ticks(tick_vals)
        cbar.set_ticklabels(tick_labels)

    cbar.ax.tick_params(which='minor', length=0)
    cbar.ax.tick_params(which='major', length=5)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200)
    else:
        plt.show()

    plt.close()



        
# ncS=Dataset(flist[0])
  
lon,lat=SST.lon.values,SST.lat.values
      
lon,lat=np.meshgrid(lon,lat)

# times=ncS['ocean_time'][:]
    
        
# phytS=ncS['phytoplankton'][0].data
# phytS[phytS>1000]=np.nan
# t=times[0]
# t_str=str(num2date(t,'seconds since 2000-1-1'))[:13]
t_str='2025-06'


sst= SST.values
sst[sst==0]=np.nan    
    
# tempS=ncS['temp'][0].data
# tempS[tempS>1000]=np.nan
draw_roms_pcolor(lon, lat, sst, timestamp=t_str, varname='Mean ERA5 SST (jun)',\
                    units=f'degree', log_scale=False, clim=(-1,32),\
                        output_path=wpth+'/'+t_str+'.png',\
                        cmap=plt.get_cmap('Spectral_r',27) )


        
             
# uS=ncS['u'][i,-22].data
# uS[uS>1000]=np.nan
# vS=ncS['v'][i,-22].data
# vS[vS>1000]=np.nan            

# speed=(uS[1:,:]**2+vS[:,1:]**2)**(1/2)


# draw_roms_pcolor(lon[1:,1:], lat[1:,1:], speed, timestamp=t_str+' (N22)', varname='Speed',\
#                  units=f'm.s', log_scale=False, clim=(0,2),\
#                      output_path=wpth+'speed/N22/'+t_str+'.png',\
#                          cmap=plt.get_cmap('Reds',27) )
            
    


       
   
   
    
    
    
    
    
    
    