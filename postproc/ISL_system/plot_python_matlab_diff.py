# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 09:54:48 2025

@author: ust21
"""

import xarray as xr
import matplotlib.pyplot as plt

# pth01="D:/shjo/ISL_TEST/test/NWP12_ini_edit2.nc"
pth01="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_zsur.nc"
pth02="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_zsur0.nc"

pth01="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_cache_zsur.nc"
pth02="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_cache_zsur0.nc"

pth01="D:/shjo/ISL_TEST/test/03_NWP12_ini_edit.nc"
pth02="D:/shjo/ISL_TEST/test/04_NWP12_ini_edit.nc"

pth01="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_test01.nc"
pth02="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_test02.nc"



pth01="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_lead.nc"
pth02="D:/shjo/ISL_TEST/test/nc/NWP12_ini_edit_cache_zsur.nc"


nc01=xr.open_dataset(pth01)['temp'][0,0].squeeze()
nc02=xr.open_dataset(pth02)['temp'][0,0].squeeze()

# nc01=nc01.where(nc02==nc02)

diff = nc01 - nc02
diff.plot(vmin=-0.02,vmax=0.02,cmap='bwr')





nc01.plot(vmin=0,vmax=3,cmap='bwr')
nc02.plot(vmin=0,vmax=3,cmap='bwr')


pth='D:/shjo/ISL_TEST/test/'





nc01=xr.open_dataset(pth01)['vbar'][0].squeeze()
nc02=xr.open_dataset(pth02)['vbar'][0].squeeze()

nc01=nc01.where(nc02==nc02)

diff = nc01 - nc02
diff.plot(vmin=-0.001,vmax=0.001,cmap='bwr')


from netCDF4 import Dataset
nc=Dataset("D:/shjo/ISL_TEST/test/NWP12_grd_edit_depth.nc")
lon,lat=nc['lon_rho'][:],nc['lat_rho'][:]

val=diff.values

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
import matplotlib.colors as mcolors  # 추가 필요

plt.rcParams.update({
    'font.size': 14,            # 기본 글꼴 크기
    'axes.titlesize': 16,       # 제목
    'axes.labelsize': 14,       # x, y축 라벨
    'xtick.labelsize': 12,      # x축 눈금
    'ytick.labelsize': 12,      # y축 눈금
    'legend.fontsize': 13,      # 범례
    'figure.titlesize': 16      # 전체 figure 제목 (suptitle)
})


fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_title(f'Temp difference surface layer', fontsize=16,pad=10)
# ax.set_extent([105, 165, 15, 50], crs=ccrs.PlateCarree())
ax.set_extent([110, 157, 15, 49], crs=ccrs.PlateCarree())

ax.coastlines(resolution='10m', linewidth=1)
ax.add_feature(cfeature.BORDERS, linestyle=':',zorder=11)
ax.add_feature(cfeature.LAND, facecolor='lightgray',zorder=10)
ax.add_feature(cfeature.OCEAN, facecolor='#ffffff')
# ax.gridlines(draw_labels=True, linestyle='--')
gl = ax.gridlines(draw_labels=True, linestyle='--',linewidth=0)
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True
# max 값 (원래 단위로)
# lonV = lon[idx][idV]
# latV = lat[idx][idV]

lonV = 0
latV = 0


# 1. log10(50)까지 색상 범위 제한
vmin = -0.02
vmax = 0.02  # ≈ 1.69897
levels = np.linspace(vmin, vmax, 33)  # 16개 구간

cmap = plt.get_cmap('bwr', len(levels) - 1)

# 2. scatter
sc = ax.pcolormesh(
    lon, lat, val,
    cmap=cmap,vmin=vmin,vmax=vmax,
   transform=ccrs.PlateCarree()
)

# 3. colorbar
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=0.8)
cbar.set_label('degree (celsius)', fontsize=12)

# ✅ tick 위치 (log10), 최대값 50까지만 표시
# tick_vals = [-2, -1.52, -1, -0.52, 0, 0.48, 1, vmax]
# tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10', '50']
# tick_vals = [-2, -1.52, -1, -0.52, 0, 0.48, 1]
# tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10']

# cbar.set_ticks(tick_vals)
# cbar.set_ticklabels(tick_labels)

# ✅ 작은 틱 제거, 큰 틱 유지
cbar.ax.tick_params(which='minor', length=0)
cbar.ax.tick_params(which='major', length=5)
plt.tight_layout()
# plt.savefig(f'D:/shjo/figs/obs_may/phyto_{t.strftime("%Y%m%d")}', dpi=200)
plt.show()



















