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

# 파일 열기
pth = 'D:/shjo/tmp/2023/'
flist=[pth+i for i in os.listdir(pth) if i.endswith('.nc')]

ref_date = np.datetime64('2000-01-01')
datasets=[]

for f in flist:
    ds = xr.open_dataset(f)
    ds = ds.expand_dims('time')

    # 날짜 추출
    fname = os.path.basename(f)
    date_str = fname.split('.')[1]  # 예: '20241231'
    dt = datetime.datetime.strptime(date_str, "%Y%m%d")
    np_time = (np.datetime64(dt) - ref_date) / np.timedelta64(1, 'D')  # float days

    ds['time'] = ('time', [np_time])
    ds['time'].attrs['units'] = 'days since 2000-01-01 00:00:00'
    datasets.append(ds)


ds_combined = xr.concat(datasets, dim='time').sel(
    lon=slice(100,170), lat=slice(55,5))

lon, lat = ds_combined.lon.values, ds_combined.lat.values

lon_m,lat_m=np.meshgrid(lon,lat)

plt.rcParams.update({
    'font.size': 14,            # 기본 글꼴 크기
    'axes.titlesize': 16,       # 제목
    'axes.labelsize': 14,       # x, y축 라벨
    'xtick.labelsize': 12,      # x축 눈금
    'ytick.labelsize': 12,      # y축 눈금
    'legend.fontsize': 13,      # 범례
    'figure.titlesize': 16      # 전체 figure 제목 (suptitle)
})



# val=ds_combined.chlor_a.loc[dict(time=slice(9283,9307))].mean(dim='time',skipna=True).values/(0.02*6.625*12)
val=ds_combined.chlor_a.mean(dim='time',skipna=True).values/(0.02*6.625*12)

# val=ds_combined.chlor_a.mean(dim='time',skipna=True).values/(0.02*6.625*12)

# val=ds_combined.chlor_a[5].values/(0.02*6.625*12)

# val2=ds_combined.chlor_a.mean(dim='time').values
# val=val2
# [115, 150, 20, 49]
fig = plt.figure(figsize=(10, 6))
# ax = plt.axes(projection=ccrs.PlateCarree())
ax = plt.axes(projection=ccrs.Mercator())

ax.set_title('Phyt (log) Observation - 0401-0430', fontsize=16,pad=10)
# ax.set_extent([105, 165, 15, 50], crs=ccrs.PlateCarree())
ax.set_extent([117, 135.5, 26, 43], crs=ccrs.PlateCarree())

ax.coastlines(resolution='10m', linewidth=1)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='#ffffff')
# ax.gridlines(draw_labels=True, linestyle='--')
gl = ax.gridlines(draw_labels=True, linestyle='--',linewidth=0)
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True
# max 값 (원래 단위로)

val_raw = val
val_clipped = np.clip(val_raw, 1e-5, None)
val_log = np.log10(val_clipped)

val_max = np.nanmax(val)
idV = np.nanargmax(val)
lonV = lon_m.flatten()[idV]
latV = lat_m.flatten()[idV]

# ax.text(
#     0.01, 0.99,
#     f'Max: {val_max:.2f} mmol/m³\n' + f'Max lon: {lonV:.2f}\n' + f'Max lat: {latV:.2f}',
#     transform=ax.transAxes,
#     fontsize=12,
#     fontweight='bold',
#     va='top',
#     ha='left',
#     bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')
# )
contour_range=[1, 16.0000, 10]
# contour_range=[-1, 1.204, 0.2]
# 1. log10(50)까지 색상 범위 제한
vmin = -1
vmax = np.log10(16)  # ≈ 1.69897
levels = np.linspace(vmin, vmax, 33)  # 16개 구간

cmap = plt.get_cmap('jet', len(levels) - 1)
cmap = plt.get_cmap('rainbow', len(levels) - 1)

norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)



# 2. scatter
sc = ax.scatter(
    lon_m.flatten(), lat_m.flatten(), c=val_log.flatten(),
    cmap=cmap, norm=norm,
    s=10, edgecolor=None, transform=ccrs.PlateCarree()
)

# 3. colorbar
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=0.8)
cbar.set_label('log(Phyt) [log(mmol m⁻³)]', fontsize=12)

# ✅ tick 위치 (log10), 최대값 50까지만 표시
# tick_vals = [-2, -1.52, -1, -0.52, 0, 0.48, 1, vmax]
# tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10', '50']
tick_vals = [-1,  -0.52, 0, 0.48, 1]
tick_labels = ['0.1',  '0.3', '1', '3', '10']

c = ax.contour(lon_m, lat_m, val,
              colors='k',
              levels=np.arange(*contour_range),
              linewidth=2.5,zorder=10, transform=ccrs.PlateCarree())
plt.clabel(c,inline=True, fontsize=15,zorder=10)

cbar.set_ticks(tick_vals)
cbar.set_ticklabels(tick_labels)

# ✅ 작은 틱 제거, 큰 틱 유지
cbar.ax.tick_params(which='minor', length=0)
cbar.ax.tick_params(which='major', length=5)
plt.tight_layout()
# plt.savefig(f'D:/shjo/ROMS_OUTS/E_test/Case_03E/figs/obs_files/phyto_{t.strftime("%Y%m%d")}', dpi=200)
plt.show()















    
    
    
    
    