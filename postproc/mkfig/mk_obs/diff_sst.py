import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.interpolate import griddata
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# ────────────────────────────────
# 사용자 설정
roms_dir = "D:/shjo/MCC/NLM_v1171_LMD/"         # ROMS 결과 폴더
ostia_file = "D:/shjo/MCC/OSTIA_SST_250601-250630.nc"        # OSTIA 파일
wpth='D:/shjo/MCC/NLM_v1171_LMD/figs/diff_sst/'
roms_varname = "temp"                # ROMS 수온 변수명
ostia_varname = "analysed_sst"                # OSTIA 수온 변수명
# ────────────────────────

# ROMS 파일 리스트 정렬
roms_files = sorted(glob.glob(os.path.join(roms_dir, "NWP15_NLM_avg_*.nc")))

# OSTIA 열기
ds_ostia = xr.open_dataset(ostia_file)
ostia_time = ds_ostia['time'].values
ostia_lat = ds_ostia['latitude'].values
ostia_lon = ds_ostia['longitude'].values
ostia_sst = ds_ostia[ostia_varname]  # 단위: Kelvin일 경우 -273.15 해줘야 함

# OSTIA 좌표 2D meshgrid
ostia_lon2d, ostia_lat2d = np.meshgrid(ostia_lon, ostia_lat)
ostia_points = np.column_stack((ostia_lon2d.ravel(), ostia_lat2d.ravel()))

# 결과 저장 리스트
bias_list = []
rmse_list = []
date_list = []


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
    ax.set_extent([105.0, 163.83, 15.0, 52.07], crs=ccrs.PlateCarree())
    # ax.set_extent([110.0, 157, 15.0, 49.04], crs=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linestyle='-.',linewidth=0)
    gl.top_labels = False
    gl.right_labels = False

    ct=ax.contour(lon2d, lat2d, val_plot,colors='k',linestyles='-',linewidths=0.8,alpha=0.85,levels=[-3,-2,-1,1,2,3])
    # ax.clabel(ct, colors='k', fontsize=12)
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


for roms_path in roms_files:
    fname = os.path.basename(roms_path)


    # ROMS 열기
    ds_roms = xr.open_dataset(roms_path)
    temp_roms = ds_roms[roms_varname].isel(s_rho=-1).squeeze().values  # 표층 temp
    lat_rho = ds_roms['lat_rho'].values
    lon_rho = ds_roms['lon_rho'].values
    mask_rho = ds_roms['mask_rho'].values  # 1: ocean, 0: land

    date_str = str(ds_roms.ocean_time.values[0])[:10]
    date_np = np.datetime64(date_str)

    # OSTIA에 해당 날짜 없으면 skip
    if date_np not in ostia_time:
        print(f"[SKIP] OSTIA 데이터 없음: {date_str}")
        continue

    # OSTIA 해당 날짜 SST 추출
    sst_ostia_day = ds_ostia[ostia_varname].sel(time=date_np).values - 273.15

    # griddata 보간 (OSTIA → ROMS grid)
    sst_interp = griddata(
        ostia_points,
        sst_ostia_day.ravel(),
        (lon_rho, lat_rho),
        method='linear'
    )

    # 유효 마스크: 바다 격자 + 유효한 SST
    valid_mask = (mask_rho == 1) & np.isfinite(sst_interp) & np.isfinite(temp_roms)

    if not np.any(valid_mask):
        print(f"[SKIP] 유효 격자 없음: {date_str}")
        continue

    # Bias / RMSE 계산
    diff = temp_roms - sst_interp
    bias = float(np.mean(diff[valid_mask]))
    rmse = float(np.sqrt(np.mean(diff[valid_mask] ** 2)))

    # 저장
    bias_list.append(bias)
    rmse_list.append(rmse)
    date_list.append(date_np)

    print(f"{date_str}  Bias={bias:.3f}°C  RMSE={rmse:.3f}°C")

    draw_roms_pcolor(lon_rho, lat_rho, diff, timestamp=date_str, varname="v1171 LMD - OSTIA",\
                        units=f'degree', log_scale=False, clim=(-3,3),\
                            output_path=wpth+date_str+'_diff_sst',\
                            cmap=plt.get_cmap('jet',27) )
    

import matplotlib.dates as mdates

# 시계열 플롯
plt.figure(figsize=(10, 4))
plt.plot(date_list, bias_list, label='Bias')
plt.plot(date_list, rmse_list, label='RMSE')
plt.title('1171 LMD vs OSTIA SST (Surface) - Daily Comparison')
plt.xlabel('Date')
plt.ylabel('°C')
plt.ylim([0.2,1.8])
# 날짜 포맷 지정 (미국식 MM/DD/YY)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d'))

# tick 간격 (2일 단위 등 조정 가능)
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(wpth+date_str+'_diff_sst_statistics')
plt.show()
