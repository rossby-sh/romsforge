import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata

# 파일 열기
obs_file = 'D:/shjo/ROMS_inputs/obs_fig/obs_sst_phyt_15km_N36.nc'
ds = nc.Dataset(obs_file)

# 변수 추출
lon = ds.variables['obs_lon'][:]
lat = ds.variables['obs_lat'][:]
val = ds.variables['obs_value'][:]
err = ds.variables['obs_error'][:]
time = ds.variables['obs_time'][:]
type = ds.variables['obs_type'][:] if 'obs_type' in ds.variables else None

# 시간 단위 변환
time_unit = ds.variables['obs_time'].units
calendar = ds.variables['obs_time'].calendar if 'calendar' in ds.variables['obs_time'].ncattrs() else 'standard'
time_dates = nc.num2date(time, time_unit, calendar=calendar)

# type == 9 (Phyt)만 선택
if type is not None:
    sst_idx = np.where(type == 9)
    lon, lat, val, err, time_dates = lon[sst_idx], lat[sst_idx], val[sst_idx], err[sst_idx], time_dates[sst_idx]

# 고정 그리드 생성
lon_min, lon_max = lon.min(), lon.max()
lat_min, lat_max = lat.min(), lat.max()
grid_x, grid_y = np.meshgrid(
    np.linspace(lon_min, lon_max, 200),
    np.linspace(lat_min, lat_max, 200)
)

# 시간 루프
unique_times = np.unique(time_dates)

for t in unique_times[90:98]:
    idx = np.where(time_dates == t)

    # log10 변환 (0 이하 값은 NaN 처리)
    val_raw = val[idx]
    val_clipped = np.clip(val_raw, 1e-5, None)
    val_log = np.log10(val_clipped)

    # 보간 (log값 기준)
    grid_val = griddata(
        (lon[idx], lat[idx]), val_log,
        (grid_x, grid_y), method='linear'
    )

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_title(f'Phyt (log10) Observation - {t.strftime("%Y-%m-%d %H:%M:%S")}', fontsize=14)
    ax.set_extent([105, 165, 15, 50], crs=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='#ffffff')
    # ax.gridlines(draw_labels=True, linestyle='--')
    gl = ax.gridlines(draw_labels=True, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    # max 값 (원래 단위로)
    val_max = np.max(val_raw)
    idV = np.argmax(val_raw)
    lonV = lon[idx][idV]
    latV = lat[idx][idV]

    ax.text(
        0.01, 0.99,
        f'Max: {val_max:.2f} mmol/m³\n' + f'Max lon: {lonV:.2f}\n' + f'Max lat: {latV:.2f}',
        transform=ax.transAxes,
        fontsize=12,
        fontweight='bold',
        va='top',
        ha='left',
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')
    )

    # scatter (log 값 기준)
    sc = ax.scatter(lon[idx], lat[idx], c=val_log, cmap='jet', s=10,
                    vmin=-2, vmax=2, edgecolor=None, transform=ccrs.PlateCarree())

    # 컬러바 (로그값 표시 + 원래 단위 텍스트 라벨)
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=.8)
    cbar.set_label('log(Phyt) [log(mmol m⁻³)]', fontsize=12)
    cbar.set_ticks([-2, -1, 0, 1,2])
    cbar.set_ticklabels(['0.01', '0.1', '1', '10','100'])

    plt.tight_layout()
    plt.show()