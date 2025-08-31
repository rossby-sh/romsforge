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
import xarray as xr

# Gpth='D:/shjo/ROMS_inputs/NWP12_grd_NWP4.nc'
# Gpth='D:/shjo/ROMS_inputs/NWP4_grd_314_10m.nc'
Gpth='D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc'

# Gpth='D:/shjo/ROMS_inputs/NWP12_grd_NWP4.nc'

Spth='D:/shjo/ROMS_OUTS/nifs_may/'

title='nifs-may'

# wpth='D:/shjo/MCC/jul025_glorys/figs/'

# flist=[Spth+i for i in os.listdir(Spth) if i.endswith('.nc')]
# flist=[Spth+i for i in os.listdir(Spth) if i.startswith('NWP12_avg_')]
flist=[Spth+i for i in os.listdir(Spth) if i.startswith('NWP15_')]

# 0226 8457

ncG=Dataset(Gpth)
lon=ncG['lon_rho'][:]
lat=ncG['lat_rho'][:]

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
    # ax.set_extent([np.min(lon2d), np.max(lon2d), np.min(lat2d), np.max(lat2d)], crs=ccrs.PlateCarree())
    ax.set_extent([105.5, 163.3, 15.0, 52.07], crs=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.',linewidth=0)
    gl.top_labels = False
    gl.right_labels = False

    ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=[4,8,12,16,20,24,28,30,32])
    ax.clabel(ct, colors='k', fontsize=12)
    # pcolormesh
    pcm = ax.pcolormesh(
        lon2d, lat2d, val_plot,
        cmap=cmap, norm=norm, shading='auto',
        transform=ccrs.PlateCarree()
    )

    # 제목
    title = f'{varname}'
    # if timestamp is not None:
    #     if hasattr(timestamp, 'strftime'):
    #         title += f' - {timestamp.strftime("%Y-%m-%d")}'
    #     else:
    #         title += f' - {timestamp}'
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




ncX=xr.open_mfdataset(flist)

temp=ncX.temp.mean(dim='ocean_time').values[-1]
phyt=ncX.phytoplankton.mean(dim='ocean_time').values[-1]


# phytS=ncS['phytoplankton'][i,-1].data



draw_roms_pcolor(lon, lat, phyt, timestamp='', varname='sst - monthly mean',\
                 units='millimole_nitrogen meter-3', log_scale=True, clim=None,\
                     output_path=None, cmap=plt.get_cmap('jet',27)) 



draw_roms_pcolor(lon, lat, temp, timestamp='', varname='sst - monthly mean',\
                    units=f'degree', log_scale=False, clim=(-1,32),\
                        output_path=None,\
                        cmap=plt.get_cmap('Spectral_r',27) )

 
   
   
   

    
    
    
    
    
    