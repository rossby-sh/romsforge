# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 12:41:49 2025

@author: ust21
"""
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
import matplotlib.colors as mcolors  # 추가 필요
import os
import xarray as xr
import datetime
import xarray as xr
import os 

pth01="D:/shjo/NWP4_NPZD/"
pth02="D:/shjo/NWP4_NPZD/NWP4_V2B/"

flist01=[pth01+i for i in os.listdir(pth01) if i.endswith(".nc")]
flist02=[pth02+i for i in os.listdir(pth02) if i.endswith(".nc")]



# ncX02.temp[0,-1].plot()
# ncX01.temp[0,-1].plot()

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
    plt.rcParams.update({
        'font.size': 14,            # 기본 글꼴 크기
        'axes.titlesize': 16,       # 제목
        'axes.labelsize': 14,       # x, y축 라벨
        'xtick.labelsize': 12,      # x축 눈금
        'ytick.labelsize': 12,      # y축 눈금
        'legend.fontsize': 13,      # 범례
        'figure.titlesize': 16      # 전체 figure 제목 (suptitle)
    })

    
    # 값 전처리
    if log_scale:
        val_plot = np.where(var2d > 1e-20, np.log10(var2d), np.nan) if log_scale else var2d.copy()
    else:
        val_plot = np.where(var2d , var2d, np.nan) if log_scale else var2d.copy()

    # val_plot = np.where(var2d , np.log10(var2d), np.nan) if log_scale else var2d.copy()

    # 색상 범위
    if clim is None:
        vmin, vmax = -3.5, 0.1
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
    ax.set_extent([110, 157, 15, 49.068], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', linewidth=1,zorder=10)
    ax.add_feature(cfeature.BORDERS, linestyle=':',zorder=10)
    ax.add_feature(cfeature.LAND, facecolor='lightgray',zorder=10)
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.',linewidth=0,zorder=10)
    gl.top_labels = False
    gl.right_labels = False

    ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=[4,8,12,16,20,24,28])
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
    # if log_scale:
    #     tick_vals = [-3, -1.52, -1, -0.52, 0.01]
    #     tick_labels = ['0.01', '0.03', '0.1', '0.3', '1']
    #     cbar.set_ticks(tick_vals)
    #     cbar.set_ticklabels(tick_labels)

    cbar.ax.tick_params(which='minor', length=0)
    cbar.ax.tick_params(which='major', length=5)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200)
    else:
        plt.show()

    plt.close()



for ii in range(12):
    
    ncX01 = xr.open_dataset(flist01[ii]).NO3[0,-1].squeeze()
    ncX02 = xr.open_dataset(flist02[ii]).NO3[0,-1].squeeze()
    
    lon01,lat01 = ncX01.lon_rho.values, ncX01.lat_rho.values
    lon02,lat02 = ncX02.lon_rho.values, ncX02.lat_rho.values
    
    
    
    val01=ncX01.values
    
    
    draw_roms_pcolor(lon01, lat01, val01, timestamp=str(ii+1)+" month", varname='NWP4 (no3)',\
                        units=f'millimole_nitrogen meter-3', log_scale=True, clim=None,\
                            output_path="D:/shjo/NWP4_NPZD/figs/no3/no3_"+str(ii+1),\
                            cmap=plt.get_cmap('rainbow',27) )
    
    
    val02=ncX02.values
    
    draw_roms_pcolor(lon02, lat02, val02, timestamp=str(ii+1)+" month", varname='NWP15 (no3)',\
                        units=f'millimole_nitrogen meter-3', log_scale=True, clim=None,\
                            output_path="D:/shjo/NWP4_NPZD/NWP4_V2B/figs/no3/no3_"+str(ii+1),\
                            cmap=plt.get_cmap('rainbow',27) )
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
