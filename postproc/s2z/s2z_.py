
#!/usr/bin/env python3
import numpy as np
import xarray as xr
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import post_utils as pu
import utils as ut   # 여기 zlevs 있음 :contentReference[oaicite:1]{index=1}


# =========================
# USER CONFIG (여기만 바꾸면 됨)
# =========================
INPUT_NC  = "/home/shjo/applications/nifs02_5km/storage/nifs02_2025_fixed_so_2/nl/nifs5km_his_9312_0005.nc"
OUTPUT_NC = "temp_z.nc"
VAR_NAME  = "temp"

ZLEVELS = np.array([-200, -150, -100, -50, -20, -10, -5, 0], dtype=float)

EXTRAP_MODE = "padding"   # "leading" or "padding"
Z_SUR = 0.0
Z_BOT = -6000.0

N_JOBS = -1
DEDUP = "jitter"          # "jitter" / "mean" / "none"

TIME_INDEX = None         # 예: 0, 전체는 None
# =========================


def _must_have(ds, names):
    for k in names:
        if k not in ds:
            raise KeyError(f"input nc에 '{k}' 없음")


def main():
    ds = xr.open_dataset(INPUT_NC)

    # ROMS 필수들 (zlevs에 필요한 것들 포함)
    need = ["h", "zeta", "hc", "Vtransform", "Vstretching", "theta_s", "theta_b"]
    _must_have(ds, need)
    if VAR_NAME not in ds:
        raise KeyError(f"input nc에 var='{VAR_NAME}' 없음")

    # time 이름 찾기
    if "ocean_time" in ds.dims:
        time_name = "ocean_time"
    elif "time" in ds.dims:
        time_name = "time"
    else:
        raise KeyError("time 차원을 찾지 못함 (ocean_time/time)")

    # 고정 필드
    h = ds["h"].values
    hc = float(ds["hc"].values)
    Vtransform = int(ds["Vtransform"].values)
    Vstretching = int(ds["Vstretching"].values)
    theta_s = float(ds["theta_s"].values)
    theta_b = float(ds["theta_b"].values)

    var = ds[VAR_NAME]

    # N 결정: 보통 s_rho 길이 = N
    if "s_rho" in ds:
        N = int(ds.dims["s_rho"])
    elif "s_rho" in var.dims:
        N = int(var.sizes["s_rho"])
    else:
        # 마지막 3D가 (time, s, y, x) 형태라는 가정이 불가하면 여기서 끊는 게 맞음
        raise KeyError("s_rho 차원을 찾지 못함 (ds['s_rho'] 또는 var.dims에 s_rho 필요)")

    # time subset
    if TIME_INDEX is None:
        zeta_all = ds["zeta"].values         # (Nt,Ny,Nx)
        var_all  = var.values                # (Nt,N,Ny,Nx) 기대
        times    = ds[time_name].values
    else:
        ti = int(TIME_INDEX)
        ds_t = ds.isel({time_name: ti})
        zeta_all = ds_t["zeta"].values[None, :, :]   # (1,Ny,Nx)
        var_all  = ds_t[VAR_NAME].values[None, ...]  # (1,N,Ny,Nx) 기대
        times    = ds_t[time_name].values[None]

    # var shape 체크
    if var_all.ndim != 4:
        raise ValueError(f"{VAR_NAME} shape 기대: (time,s_rho,eta,xi), 실제: {var_all.shape}")

    Nt, Ns, Ny, Nx = var_all.shape
    if Ns != N:
        raise ValueError(f"N mismatch: N={N}, var vertical={Ns}")

    Nz = ZLEVELS.size
    out = np.empty((Nt, Nz, Ny, Nx), dtype=float)

    # ---- 핵심: time마다 utils.zlevs로 z_rho 생성 :contentReference[oaicite:2]{index=2}
    # utils.zlevs 시그니처:
    # zlevs(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, zeta)
    # igrid=1 => rho point
    for t in range(Nt):
        zeta2d = zeta_all[t, :, :]
        z_rho = ut.zlevs(
            Vtransform, Vstretching, theta_s, theta_b,
            hc, N, 1, h, zeta2d
        )  # (N,Ny,Nx) 로 나오는게 정상 :contentReference[oaicite:3]{index=3}

        v3d = var_all[t, :, :, :]  # (N,Ny,Nx)

        if EXTRAP_MODE == "padding":
            out[t] = pu.interpolate_s_to_zlevels(
                z_rho, v3d, ZLEVELS,
                n_jobs=N_JOBS,
                dedup=DEDUP,
                extrap_mode="padding",
                zsur=Z_SUR,
                zbot=Z_BOT
            )
        elif EXTRAP_MODE == "leading":
            out[t] = pu.interpolate_s_to_zlevels(
                z_rho, v3d, ZLEVELS,
                n_jobs=N_JOBS,
                dedup=DEDUP,
                extrap_mode="leading"
            )
        else:
            raise ValueError("EXTRAP_MODE는 'leading' 또는 'padding'만 가능")

    # dims 이름은 var의 마지막 2개를 그대로 씀
    eta_name = var.dims[-2]
    xi_name  = var.dims[-1]

    ds_out = xr.Dataset()
    ds_out[time_name] = xr.DataArray(times, dims=(time_name,))
    ds_out["z"] = xr.DataArray(ZLEVELS, dims=("z",), attrs={"units": "m", "positive": "up"})
    ds_out[VAR_NAME] = xr.DataArray(
        out,
        dims=(time_name, "z", eta_name, xi_name),
        attrs=dict(var.attrs) if hasattr(var, "attrs") else {}
    )

    # 있으면 같이 복사
    for name in ["lon_rho", "lat_rho", "mask_rho"]:
        if name in ds:
            ds_out[name] = ds[name]

    ds_out["h"] = xr.DataArray(h, dims=(eta_name, xi_name), attrs={"units": "m"})
    ds_out.to_netcdf(OUTPUT_NC)

    ds.close()
    print(f"written: {OUTPUT_NC}")


if __name__ == "__main__":
    main()
