# --- [00] Imports and path setup ---
import sys
import os
import glob
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl
import time 
from pathlib import Path

# --- [01] Load configuration and input metadata ---
cfg  = tl.parse_config("./config_proc_nft.yaml")
grd  = tl.load_roms_grid(cfg.grdname)

base_dir = Path(getattr(cfg, "base_dir", ".")).expanduser().resolve()
case = getattr(cfg, "case", "case")

outdir = base_dir / case / "figures"
datadir = base_dir / case 
outdir.mkdir(parents=True, exist_ok=True)
for sub in ["temp", "salt", "phyt", "chl"]:
    (outdir / sub).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors


title = case

chl_bool = 1

# flist=[Spth+i for i in os.listdir(Spth) if i.endswith('.nc')]
flist=[str(datadir)+'/'+i for i in os.listdir(str(datadir)) if i.endswith('.nc')]


ncG=Dataset(cfg.grdname)
lon=ncG['lon_rho'][:]
lat=ncG['lat_rho'][:]

def draw_roms_pcolor(lon2d, lat2d, var2d, timestamp=None,
                     varname='Variable', units='', log_scale=False,
                     clim=None, output_path=None, cmap='jet',levels=[]):
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
        # val_plot = np.where(var2d > 1e-30, np.log10(var2d), np.nan) if log_scale else var2d.copy()
        eps = 1e-20  # 너무 작은 값 방지
        val_plot = np.log10(np.clip(var2d, eps, None))
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
    ax.set_extent([np.min(lon2d), np.max(lon2d), np.min(lat2d), np.max(lat2d)], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.',linewidth=0)
    gl.top_labels = False
    gl.right_labels = False

    # ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=[4,8,12,16,20,24,28,30,32])
    ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=levels)

    ax.clabel(ct, colors='k', fontsize=12)
    # pcolormesh
    # pcm = ax.contourf(
    #     lon2d, lat2d, val_plot,
    #     cmap=cmap, norm=norm, shading='auto',extend="both",levels=[-2, -1.52, -1, -0.52, 0, 0.48, 1],
    #     transform=ccrs.PlateCarree()
    # )
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
    label = f'[{units}]'
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


for f in flist:

    with Dataset(f) as ncS:
        
        times=ncS['ocean_time'][:]
           
        for i in range(len(times)):
            t=times[i]
            t_str=str(num2date(t,'seconds since 2000-1-1'))[:10]
            print(t_str)
            
            phytS=ncS['phytoplankton'][i].data
            phytS[phytS>1000]=np.nan
            draw_roms_pcolor(lon, lat, phytS, timestamp=t_str, varname='phyt '+title,\
                             units='millimole_nitrogen meter-3', log_scale=True, clim=None,\
                                 output_path=str(outdir)+'/phyt/'+t_str+'.png', cmap=plt.get_cmap('jet',27))
            

            tempS=ncS['temp'][i].data
            tempS[tempS>1000]=np.nan
            draw_roms_pcolor(lon, lat, tempS, timestamp=t_str, varname='SST '+title,\
                                units=f'degree', log_scale=False, clim=(-1,32),\
                                    output_path=str(outdir)+'/temp/'+t_str+'.png',\
                                    cmap=plt.get_cmap('Spectral_r',27),levels=[4,8,12,16,20,24,28,30,32])
            
            if chl_bool:
                chloS=ncS['chlorophyll'][i].data
                chloS[chloS>1000]=np.nan   
            else:
                chloS=phytS*(0.02*6.625*12)
            draw_roms_pcolor(lon, lat, chloS, timestamp=t_str, varname='chl '+title,\
                             units='millimole_nitrogen meter-3', log_scale=True, clim=None,\
                                 output_path=str(outdir)+'/chl/'+t_str+'.png', cmap=plt.get_cmap('jet',27))
            
            
            # if np.nanmin(chloS)<0:
            #     print(np.nanmin(chloS))
            
            
            saltS=ncS['salt'][i].data
            saltS[saltS>1000]=np.nan
            draw_roms_pcolor(lon, lat, saltS, timestamp=t_str, varname='SSS '+title,\
                                units=f'psu', log_scale=False, clim=(28,35.5),\
                                    output_path=str(outdir)+'/salt/'+t_str+'.png',\
                                    cmap=plt.get_cmap("jet"),levels=[28,30,32,33,34,34.5,35]) 
            

            # uS=ncS['u'][i,-22].data
            # uS[uS>1000]=np.nan
            # vS=ncS['v'][i,-22].data
            # vS[vS>1000]=np.nan            
            
            # speed=(uS[1:,:]**2+vS[:,1:]**2)**(1/2)

            # draw_roms_pcolor(lon[1:,1:], lat[1:,1:], speed, timestamp=t_str+' (N22)', varname='Speed',\
            #                  units=f'm.s', log_scale=False, clim=(0,2),\
            #                      output_path=outdir+'speed/N22/'+t_str+'.png',\
            #                          cmap=plt.get_cmap('Reds',27) )
                        
    # cmocean.cm.haline

    
            TIME = [f"{h:02d}:00" for h in np.arange(0,24,3)]

    
    
    
    
    
    
    


