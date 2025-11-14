
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_zeta_uv.py  zeta(해수면 높이) + UV 벡터 표시
  - zeta: cmap='balance', vmin=-0.5, vmax=2
  - UV: 정규화 방향 화살표, 색은 속도(rainbow)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset, num2date
import cmocean

# ------------------ 내부 보조 함수 ------------------
def uv_to_rho_interior(u_u, v_v, zeta_rho, lon_rho, lat_rho):
    """
    C-그리드(u,v) → rho 내부로 변환
    u_u: (eta, xi_u)
    v_v: (eta_v, xi)
    반환: (zeta_i, lon_i, lat_i, u_i, v_i)  모두 동일 shape (eta-2, xi-2)
    """
    u_c = 0.5 * (u_u[:, :-1] + u_u[:, 1:])   # xi 방향 평균
    v_c = 0.5 * (v_v[:-1, :] + v_v[1:, :])   # eta 방향 평균

    z_i   = zeta_rho[1:-1, 1:-1]
    lon_i = lon_rho [1:-1, 1:-1]
    lat_i = lat_rho [1:-1, 1:-1]
    u_i   = u_c[1:-1, :]
    v_i   = v_c[:, 1:-1]

    eta_min = min(z_i.shape[0], u_i.shape[0], v_i.shape[0])
    xi_min  = min(z_i.shape[1], u_i.shape[1], v_i.shape[1])
    return (
        z_i[:eta_min, :xi_min],
        lon_i[:eta_min, :xi_min],
        lat_i[:eta_min, :xi_min],
        u_i[:eta_min, :xi_min],
        v_i[:eta_min, :xi_min],
    )

def draw_zeta_uv(lon2d, lat2d, zeta2d, u_u, v_v,
                 timestamp=None, output_path=None,
                 clim=(-0.5, 2.0), stride=6):
    """zeta(pcolormesh) + uv(quiver)"""
    zI, lonI, latI, uI, vI = uv_to_rho_interior(u_u, v_v, zeta2d, lon2d, lat2d)

    speed = np.sqrt(uI**2 + vI**2)
    mag = np.maximum(speed, 1e-12)
    uN, vN = uI / mag, vI / mag

    cmap_zeta = cmocean.cm.balance
    cmap_uv   = plt.get_cmap("rainbow")

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lonI.min(), lonI.max(), latI.min(), latI.max()], crs=ccrs.PlateCarree())
    ax.coastlines(resolution="10m", linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.6)
    ax.add_feature(cfeature.LAND, facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN, facecolor="white")
    gl = ax.gridlines(draw_labels=True, linestyle="-.", linewidth=0)
    gl.top_labels = False
    gl.right_labels = False

    pcm = ax.pcolormesh(lonI, latI, zI, cmap=cmap_zeta,
                        shading="auto", vmin=clim[0], vmax=clim[1],
                        transform=ccrs.PlateCarree())

    lonQ, latQ = lonI[::stride, ::stride], latI[::stride, ::stride]
    uQ, vQ, sQ = uN[::stride, ::stride], vN[::stride, ::stride], speed[::stride, ::stride]

    mask = np.isfinite(lonQ) & np.isfinite(latQ) & np.isfinite(uQ) & np.isfinite(vQ) & np.isfinite(sQ)
    lonQ, latQ, uQ, vQ, sQ = lonQ[mask], latQ[mask], uQ[mask], vQ[mask], sQ[mask]

    Q = ax.quiver(lonQ, latQ, uQ, vQ, sQ, cmap=cmap_uv, scale=40, width=0.0025,
                  transform=ccrs.PlateCarree(), pivot="middle")

    ax.set_title(f"Zeta + Surface Velocity - {timestamp}" if timestamp else "Zeta + Surface Velocity", pad=10)
    cbar1 = plt.colorbar(pcm, ax=ax, pad=0.02, shrink=0.82)
    cbar1.set_label("zeta [m]")
    cbar2 = plt.colorbar(Q, ax=ax, pad=0.08, shrink=0.82)
    cbar2.set_label("speed [m/s]")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=200)
    else:
        plt.show()
    plt.close()


# ------------------ 실행부 ------------------
if __name__ == "__main__":
    # 예시: 직접 경로 수정
    ncpath = "/home/shjo/warehouse/nifs01/leev2_tak_da8ol4_ssh_sst_kodc/surface_monthly_avg_leev2_tak_da8ol4_ssh_sst_kodc.nc"
    with Dataset(ncpath) as nc:
        lon = nc["lon_rho"][:]
        lat = nc["lat_rho"][:]
        times = nc["ocean_time"][:]

        for i in range(len(times)):
            t = num2date(times[i], "seconds since 2000-01-01 00:00:00")
            t_str = str(t)[:10]
            print(f"[{t_str}] plotting {os.path.basename(ncpath)}")

            zeta = nc["zeta"][i, :, :].data
            u2d = nc["u"][i, :, :].data
            v2d = nc["v"][i, :, :].data

            out_path = f"{t_str}_zeta_uv.png"
            draw_zeta_uv(lon, lat, zeta, u2d, v2d, timestamp=t_str, output_path=out_path)
