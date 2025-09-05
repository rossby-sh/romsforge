import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
import matplotlib.colors as mcolors  # 추가 필요

# 파일 열기
obs_file = 'D:/shjo/ROMS_OUTS/E_test/roms_obs_phyt_30km_N36_250501_250531_30Elimited002_.nc'
ds = nc.Dataset(obs_file)

plt.rcParams.update({
    'font.size': 14,            # 기본 글꼴 크기
    'axes.titlesize': 16,       # 제목
    'axes.labelsize': 14,       # x, y축 라벨
    'xtick.labelsize': 12,      # x축 눈금
    'ytick.labelsize': 12,      # y축 눈금
    'legend.fontsize': 13,      # 범례
    'figure.titlesize': 16      # 전체 figure 제목 (suptitle)
})

# 변수 추출
lon = ds.variables['obs_lon'][:]
lat = ds.variables['obs_lat'][:]
val = ds.variables['obs_value'][:]
err = ds.variables['obs_error'][:]
time = ds.variables['obs_time'][:]
type = ds.variables['obs_type'][:] if 'obs_type' in ds.variables else None

# 시간 단위 변환
time_unit = ds.variables['obs_time'].units
time_unit = 'days since 2000-1-1'
calendar = ds.variables['obs_time'].calendar if 'calendar' in ds.variables['obs_time'].ncattrs() else 'standard'
time_dates = nc.num2date(time, time_unit, calendar=calendar)

# type == 9 (Phyt)만 선택
if type is not None:
    sst_idx = np.where(type == 9)
    lon, lat, val, err, time_dates = lon[sst_idx], lat[sst_idx], val[sst_idx], err[sst_idx], time_dates[sst_idx]


idS=np.where(lat<=25)

    
lonS, latS, valS, errS, time_datesS=lon[idS], lat[idS], val[idS], err[idS], time_dates[idS]
    
idE=np.where(lonS>=127)
lonS, latS, valS, errS, time_datesS=lonS[idE], latS[idE], valS[idE], errS[idE], time_datesS[idE]


import pandas as pd



# df = pd.DataFrame(columns=['mean', 'min', 'max','mean(log)', 'min(log)', 'max(log)',], index=pd.date_range('2025-05-01', '2025-05-31'))
df = pd.DataFrame(columns=['mean', 'min', 'max','mean(log)', 'min(log)', 'max(log)',], index=range(31))


# 시간 루프
unique_times = np.unique(time_dates)

for t,n in zip(unique_times,range(31)):
    idx = np.where(time_datesS == t)

    # log10 변환 (0 이하 값은 NaN 처리)
    val_raw = valS[idx]
    val_clipped = np.clip(val_raw, 1e-5, None)
    val_log = np.log10(val_clipped)
    print(np.min(errS[idx]))


    MEAM=np.mean(val_clipped)
    MAX=np.max(val_clipped)
    MIN=np.min(val_clipped)

    MEAM_L=np.mean(val_log)
    MAX_L=np.max(val_log)
    MIN_L=np.min(val_log)
    
    df.loc[n] = [MEAM, MIN, MAX, MEAM_L, MIN_L, MAX_L]

df.index=pd.date_range('2025-05-01', '2025-05-31')

# df.to_csv('D:/shjo/ROMS_OUTS/E_test/Case_03E/figs/obs_values/25S_127E.csv')

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_title(f'Phyt (log) Observation - {t.strftime("%Y-%m-%d")}', fontsize=16,pad=10)
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

    # 1. log10(50)까지 색상 범위 제한
    vmin = -2
    vmax = np.log10(20)  # ≈ 1.69897
    levels = np.linspace(vmin, vmax, 33)  # 16개 구간

    cmap = plt.get_cmap('jet', len(levels) - 1)
    norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)

    # 2. scatter
    sc = ax.scatter(
        lonS[idx], latS[idx], c=val_log,
        cmap=cmap, norm=norm,
        s=10, edgecolor=None, transform=ccrs.PlateCarree()
    )

    # 3. colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=0.8)
    cbar.set_label('log(Phyt) [log(mmol m⁻³)]', fontsize=12)

    # ✅ tick 위치 (log10), 최대값 50까지만 표시
    # tick_vals = [-2, -1.52, -1, -0.52, 0, 0.48, 1, vmax]
    # tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10', '50']
    tick_vals = [-2, -1.52, -1, -0.52, 0, 0.48, 1]
    tick_labels = ['0.01', '0.03', '0.1', '0.3', '1', '3', '10']
    
    cbar.set_ticks(tick_vals)
    cbar.set_ticklabels(tick_labels)

    # ✅ 작은 틱 제거, 큰 틱 유지
    cbar.ax.tick_params(which='minor', length=0)
    cbar.ax.tick_params(which='major', length=5)
    plt.tight_layout()
    # plt.savefig(f'D:/shjo/ROMS_OUTS/E_test/Case_03E/figs/obs_files/phyto_{t.strftime("%Y%m%d")}', dpi=200)
    plt.show()
    
    













    
    
    
    
    