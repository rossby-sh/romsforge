# -*- coding: utf-8 -*-
"""
Created on Wed May 21 08:29:33 2025

@author: ust21
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import numpy as np

from netCDF4 import Dataset
import cmocean

Gpth='D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc'
Spth='D:/shjo/stdfile_plot/NIFS_NWP_Ana_202504.nc'



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
    val_plot = np.where(var2d > 1e-5, np.log10(var2d), np.nan) if log_scale else var2d.copy()

    # 색상 범위
    if clim is None:
        vmin, vmax = np.nanmin(val_plot), np.nanmax(val_plot)
    else:
        vmin, vmax = clim

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
    ax.set_extent([np.min(lon2d), np.max(lon2d), np.min(lat2d), np.max(lat2d)], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.')
    gl.top_labels = False
    gl.right_labels = False

    # pcolormesh
    pcm = ax.pcolormesh(
        lon2d, lat2d, val_plot,
        cmap=cmap, norm=norm, shading='auto',
        transform=ccrs.PlateCarree()
    )

    # 제목
    title = f'{varname} (Model)'
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
        tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10', '50']
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





ncS=Dataset(Spth)
lon=ncS['lon'][:]
lat=ncS['lat'][:]
tempS=ncS['temp'][0]
saltS=ncS['salt'][0]
zetaS=ncS['zeta'][:]

tempS=ncS['temp'][0]


NO3S=ncS['NO3'][0]
phytS=ncS['phytoplankton'][0]
zoopS=ncS['zooplankton'][0]


# Temperature
draw_roms_pcolor(lon, lat, tempS, timestamp=None, varname='Surface temperature',\
                 units='celcius', log_scale=False, clim=[0,25], output_path=None, cmap='RdYlBu_r')
    
# saltS
draw_roms_pcolor(lon, lat, saltS, timestamp=None, varname='Surface salinity',\
                 units='psu', log_scale=False, clim=[32,35], output_path=None, cmap='jet')
    
# zeta
draw_roms_pcolor(lon, lat, zetaS, timestamp=None, varname='zeta',\
                 units='m', log_scale=False, clim=[-1,1], output_path=None, cmap=cmocean.cm.balance) 
    
# NO3S
draw_roms_pcolor(lon, lat, NO3S, timestamp=None, varname='Surface NO3',\
                 units='millimole_nitrogen meter-3', log_scale=False, clim=None, output_path=None, cmap='jet') 
    
# phytS
draw_roms_pcolor(lon, lat, phytS, timestamp=None, varname='Surface phyt',\
                 units='millimole_nitrogen meter-3', log_scale=True, clim=None, output_path=None, cmap='jet') 
    
draw_roms_pcolor(lon, lat, phytS, timestamp=None, varname='Surface phyt',\
                 units='millimole_nitrogen meter-3', log_scale=False, clim=None, output_path=None, cmap='jet') 
    
# zoopS
draw_roms_pcolor(lon, lat, zoopS, timestamp=None, varname='Surface zoop',\
                 units='millimole_nitrogen meter-3', log_scale=False, clim=None, output_path=None, cmap='jet') 
    
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    