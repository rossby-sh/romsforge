# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:24:57 2025

@author: ust21
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
from netCDF4 import Dataset
import cmocean

npth00 = "D:/shjo/LTE/NWP12_grd_NWP4.nc"
npth01 = "D:/shjo/LTE/NWP12_grd_edit_depth.nc"
npth02 = "D:/shjo/LTE/NWP4_grd_3_10m_LP_subsampled.nc"
npth03 = "D:/shjo/LTE/NWP4_grd_3_10m_LP.nc"

ncG00=Dataset(npth00)
ncG01=Dataset(npth01)
ncG02=Dataset(npth02,"r")
ncG03=Dataset(npth03,"r")

lon01, lat01 = ncG01['lon_rho'][:].data, ncG01['lat_rho'][:].data 
lon02, lat02 = ncG02['lon_rho'][:].data, ncG02['lat_rho'][:].data 

mask_rho01=ncG01['mask_rho'][:]
mask_rho02=ncG02['mask_rho'][:]

h01 = ncG01["h"][:].data

h_sub = h01[::3,::3]
h_edt = h01[:]
h_ori=ncG00['h'][:]
h_4=ncG03['h'][:]

# ncG02["h"][:] = h_sub

ncG02.close()
ncG01.close()
ncG00.close()
ncG03.close()


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
    # ax.set_extent([105.0, 163.83, 15.0, 52.07], crs=ccrs.PlateCarree())
    # ax.set_extent([115.0, 150.0, 20.0, 49.0], crs=ccrs.PlateCarree())
    ax.set_extent([145.0, 157.0, 20.0, 30.0], crs=ccrs.PlateCarree())
    ax.set_extent([np.min(lon2d), np.max(lon2d), np.min(lat2d), np.max(lat2d)], crs=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.',linewidth=0)
    gl.top_labels = False
    gl.right_labels = False

    # ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=[4,8,12,16,20,24,25,26,27,28])
    # ax.clabel(ct, colors='k', fontsize=12)
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



draw_roms_pcolor(lon02, lat02, h_sub, timestamp='', varname='subsampled',\
                    units=f'm', log_scale=False, clim=(7000,0),\
                        output_path=None,\
                        cmap=cmocean.cm.deep )

draw_roms_pcolor(lon02, lat02, h_4, timestamp='', varname='NWP4',\
                    units=f'm', log_scale=False, clim=(7000,0),\
                        output_path=None,\
                        cmap=cmocean.cm.deep )


draw_roms_pcolor(lon01, lat01, h, timestamp='', varname='NWP12_modi',\
                    units=f'm', log_scale=False, clim=(7000,0),\
                        output_path=None,\
                        cmap=cmocean.cm.deep )
    
draw_roms_pcolor(lon01, lat01, h_ori, timestamp='', varname='NWP12',\
                    units=f'm', log_scale=False, clim=(7000,0),\
                        output_path=None,\
                        cmap=cmocean.cm.deep )





import numpy as np

def compute_rfactor(h: np.ndarray, mask: np.ndarray | None = None):
    """
    ROMS-style Topography R-factor 계산
    h : 2D numpy array (positive down depth)
    mask : 2D boolean array (True=ocean), optional

    return:
        R, Rx, Ry (모두 2D array, shape=h.shape)
    """
    if h.ndim != 2:
        raise ValueError("h는 반드시 2D array여야 함")

    finite = np.isfinite(h)
    pos = h > 0
    valid = finite & pos if mask is None else (mask & finite & pos)
    # valid = h

    h = np.where(valid, h, np.nan)

    # x방향 경사
    dx = np.abs(h[:, 1:] - h[:, :-1])
    sx = h[:, 1:] + h[:, :-1]
    Rx = dx / sx

    # y방향 경사
    dy = np.abs(h[1:, :] - h[:-1, :])
    sy = h[1:, :] + h[:-1, :]
    Ry = dy / sy

    # cell-centered로 확장
    Rx_full = np.full_like(h, np.nan)
    Ry_full = np.full_like(h, np.nan)

    Rx_full[:, 1:] = Rx
    Rx_full[:, :-1] = np.fmax(Rx_full[:, :-1], Rx)

    Ry_full[1:, :] = Ry
    Ry_full[:-1, :] = np.fmax(Ry_full[:-1, :], Ry)

    R = np.fmax(Rx_full, Ry_full)

    return R, Rx_full, Ry_full

# mask = (h != ).astype(int)   # 바다는 1, 육지는 0
# mask_bool = (h != 5)   # 바다는 True, 육지는 False




mask_bool01 = mask_rho01==True
mask_bool02 = mask_rho02==True

R_edt,_,_=compute_rfactor(h_edt,mask_bool01)
R_sub,_,_=compute_rfactor(h_sub,mask_bool02)
R_4,_,_=compute_rfactor(h_4,mask_bool02)


# a[a>0.15]=np.nan
# plt.pcolor(R_ori,vmin=0,vmax=0.2); plt.colorbar()





import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def plot_rfactor_easy(R, mask=None, invert_y=True):
    """
    R: 2D float array (R-factor)
    mask: 2D bool (True=바다, False=육지) 있으면 육지는 회색 처리
    """
    R = np.array(R, dtype=float)

    # 육지는 NaN으로 가려두기 (있으면)
    if mask is not None:
        R = np.where(mask, R, np.nan)

    # 등급 숫자화: 0(안전) / 1(주의) / 2(위험)
    Z = np.full(R.shape, np.nan, dtype=float)
    Z[R <= 0.2] = 0.0
    Z[(R > 0.2) & (R <= 0.3)] = 1.0
    Z[R > 0.3] = 2.0

    # Z[R <= 0.1] = 0.0
    # Z[(R > 0.1) & (R <= 0.2)] = 1.0
    # Z[R > 0.2] = 2.0

    # 컬러맵: 0→파랑, 1→주황, 2→빨강, NaN→회색
    cmap = ListedColormap(["blue", "orange", "red"])
    cmap.set_bad("lightgray")

    # imshow로 간단히: vmin=-0.5, vmax=2.5면 0/1/2가 딱 떨어짐
    plt.figure(figsize=(7,6))
    im = plt.imshow(Z, cmap=cmap, vmin=-0.5, vmax=3, origin="upper")
    if invert_y:
        plt.gca().invert_yaxis()  # 그리드 위쪽이 y=0이면 보기 좋게 뒤집기

    # 간단한 범례
    from matplotlib.patches import Patch
    legend = [
        Patch(color="blue", label="≤ 0.2 (stable)"),
        Patch(color="orange", label="0.2–0.3 (caution)"),
        Patch(color="red", label="> 0.3 (risky)")
    ]
    plt.legend(handles=legend, loc="upper left", frameon=True)

    plt.title("R-factor Stability")
    plt.xlabel("X index"); plt.ylabel("Y index")
    plt.tight_layout(); plt.show()

# 예시
if __name__ == "__main__":
    plot_rfactor_easy(R_sub)  # mask_rho 같은 게 있으면 plot_rfactor_easy(R, mask_rho)

    
plot_rfactor_easy(R_edt)  # mask_rho 같은 게 있으면 plot_rfactor_easy(R, mask_rho)
plot_rfactor_easy(R_4)  # mask_rho 같은 게 있으면 plot_rfactor_easy(R, mask_rho)
    
np.nanmax(R_sub)
    

    
    
    
    
    
    
    
    
    
    
    
    







