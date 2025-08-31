# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 09:48:04 2025

@author: ust21
"""

import xarray as xr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
import matplotlib.colors as mcolors  # 추가 필요
from matplotlib.ticker import FixedLocator  # ✅ 여기서 가져오기
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os

nrpth='C:/Users/ust21/shjo/MCC/X_nifs_ini_20250701.nc'
ncA=xr.open_dataset(nrpth)
ncG=xr.open_dataset('D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc')

# rpth='C:/Users/ust21/shjo/LTE/fwd/'
# POST_list=[rpth+i for i in os.listdir(rpth) if i.endswith("001.nc")]
# PRIOR_list=[rpth+i for i in os.listdir(rpth) if i.endswith("000.nc")]

# POST01=xr.open_mfdataset(POST_list[0]).salt.loc[dict(ocean_time=slice("2023-04-01 00","2023-04-07 00"))]
# POST02=xr.open_mfdataset(POST_list[1]).salt.loc[dict(ocean_time=slice("2023-04-07 02","2023-04-15 00"))]
# POST03=xr.open_mfdataset(POST_list[2]).salt.loc[dict(ocean_time=slice("2023-04-15 02","2023-04-23 00"))]
# POST04=xr.open_mfdataset(POST_list[3]).salt.loc[dict(ocean_time=slice("2023-04-23 02","2023-04-30 23"))]

# POST=xr.concat([POST01,POST02,POST03,POST04],dim='ocean_time').mean(dim='ocean_time')



lon_rho,lat_rho=ncG.lon_rho.values,ncG.lat_rho.values


val=ncA['zooplankton'].values[0,-1]

lon_ticks = np.arange(115, 157, 10)  # 117, 122, 127, 132
lat_ticks = np.arange(15, 46, 5)    # 26, 31, 36, 41

plt.rcParams.update({
    'font.weight': 'bold',      # 기본 폰트를 bold로
    'axes.labelweight': 'bold', # x/y축 레이블
    'axes.titleweight': 'bold', # 제목
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})

# =============================================================================
# plot post
# =============================================================================
contour_range=np.arange(28,36,2)
contour_range2=np.arange(28,36,2)

# contour_range=[-1, 1.204, 0.2]
# 1. log10(50)까지 색상 범위 제한
vmin = 0
vmax = 1  # ≈ 1.69897
levels = np.linspace(vmin, vmax, 64)  # 16개 구간
cmap = plt.get_cmap('rainbow', len(levels) - 1)

val_plot = np.where(val > 1e-20, np.log10(val), np.nan) if True else val.copy()
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# figure & axes
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.Mercator())
ax.set_extent([105, 163, 15, 50], crs=ccrs.PlateCarree())
ax.set_title('zooplankton', fontsize=16, pad=10)

# 1. OCEAN / LAND 배경
ax.add_feature(cfeature.OCEAN, facecolor='#ffffff', zorder=0)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=10)

# 2. 데이터 contourf
sc = ax.contourf(
    lon_rho, lat_rho, val_plot,
    vmin=vmin, vmax=vmax, levels=levels, extend='both',
    cmap=cmap, transform=ccrs.PlateCarree(),
    zorder=2, norm=norm
)

# 3. 경계선과 해안선 덮어쓰기
ax.coastlines(resolution='10m', linewidth=1, color='k', zorder=13)
ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='k', zorder=13)

# 4. contour 선과 레이블
# c = ax.contour(
#     lon_rho, lat_rho, val_plot,
#     colors='k',
#     levels=contour_range,
#     linewidths=1.3,
#     zorder=4, transform=ccrs.PlateCarree()
# )
# plt.clabel(c, inline=True, fontsize=10, zorder=5)

# 5. colorbar
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=1)
cbar.set_label('', fontsize=12)
# cbar.set_ticks(contour_range2)
# cbar.set_ticklabels(contour_range2)
cbar.ax.tick_params(which='minor', length=0)
cbar.ax.tick_params(which='major', length=5)

tick_vals = [-1.52, -1, -0.52, 0, 0.48, 1]
tick_labels = ['0.03', '0.1', '0.3', '1', '3', '10']
cbar.set_ticks(tick_vals)
cbar.set_ticklabels(tick_labels)

# 6. 그리드와 눈금 (117, 26도 표시되도록)
gl = ax.gridlines(draw_labels=False, linestyle='--', color='gray', alpha=0.7,linewidth=0)
# draw_labels=False → 눈금 포맷을 직접 지정

ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())

plt.tight_layout()
plt.show()













