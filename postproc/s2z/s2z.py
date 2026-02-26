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
OUTPUT_NC = "temp_z2.nc"

VARS_RHO = ["temp", "salt"]     # (time, s_rho, eta_rho, xi_rho)
VAR_U    = "u"                  # (time, s_rho, eta_u,  xi_u)
VAR_V    = "v"                  # (time, s_rho, eta_v,  xi_v)

ZLEVELS = np.array([-200, -150, -100, -50, -20, -10, -5, 0], dtype=float)

EXTRAP_MODE = "padding"   # "leading" or "padding"
Z_SUR = 0.0
Z_BOT = -6000.0

N_JOBS = -1
DEDUP = "none"          # "jitter" / "mean" / "none"

TIME_INDEX = None         # 예: 0, 전체는 None
# =========================

def _must_have(ds, names):
    for k in names:
        if k not in ds:
            raise KeyError(f"input nc에 '{k}' 없음")


def _coord_sanitized(da: xr.DataArray) -> xr.DataArray:
    """
    xarray merge 꼬임 방지용:
    da가 원래 갖고 있던 coords 메타를 버리고 values+dims만 살려서 새로 만든다.
    """
    return xr.DataArray(da.values, dims=da.dims, attrs=dict(da.attrs))


def _infer_time_name(ds):
    if "ocean_time" in ds.dims:
        return "ocean_time"
    if "time" in ds.dims:
        return "time"
    raise KeyError("time 차원을 찾지 못함 (ocean_time/time)")


def _interp_one_time(z3d, v3d):
    if EXTRAP_MODE == "padding":
        return pu.interpolate_s_to_zlevels(
            z3d, v3d, ZLEVELS,
            n_jobs=N_JOBS,
            dedup=DEDUP,
            extrap_mode="padding",
            zsur=Z_SUR,
            zbot=Z_BOT
        )
    elif EXTRAP_MODE == "leading":
        return pu.interpolate_s_to_zlevels(
            z3d, v3d, ZLEVELS,
            n_jobs=N_JOBS,
            dedup=DEDUP,
            extrap_mode="leading"
        )
    else:
        raise ValueError("EXTRAP_MODE는 'leading' 또는 'padding'만 가능")


def _make_z(ds, h_eta_xi, zeta_eta_xi, N, igrid):
    """
    utils.zlevs()는 (xi, eta) 순서를 가정하고 평균(igrid=3/4)을 내므로,
    입력을 transpose해서 넣고, 출력도 다시 (eta, xi)로 되돌린다.
    """
    Vtransform  = int(ds["Vtransform"].values)
    Vstretching = int(ds["Vstretching"].values)
    theta_s = float(ds["theta_s"].values)
    theta_b = float(ds["theta_b"].values)
    hc = float(ds["hc"].values)

    # (eta, xi) -> (xi, eta)
    z_xi_eta = ut.zlevs(
        Vtransform, Vstretching, theta_s, theta_b,
        hc, N, igrid,
        h_eta_xi.T, zeta_eta_xi.T
    )  # (N, xi, eta) 형태로 나옴

    # (N, xi, eta) -> (N, eta, xi)
    z_eta_xi = np.transpose(z_xi_eta, (0, 2, 1))
    return z_eta_xi

def main():
    ds = xr.open_dataset(INPUT_NC)

    # zlevs에 필요한 기본들
    need = ["h", "zeta", "hc", "Vtransform", "Vstretching", "theta_s", "theta_b"]
    _must_have(ds, need)

    time_name = _infer_time_name(ds)

    # N (FutureWarning 회피: ds.dims 대신 ds.sizes)
    if "s_rho" in ds.sizes:
        N = int(ds.sizes["s_rho"])
    else:
        raise KeyError("s_rho dimension이 없음")

    h = ds["h"].values

    # time subset
    if TIME_INDEX is None:
        zeta_all = ds["zeta"].values                     # (Nt,eta_rho,xi_rho)
        times = ds[time_name].values
        t_indices = range(zeta_all.shape[0])
    else:
        ti = int(TIME_INDEX)
        zeta_all = ds["zeta"].isel({time_name: ti}).values[None, :, :]
        times = ds[time_name].isel({time_name: ti}).values[None]
        t_indices = range(1)

    Nt = len(t_indices)
    Nz = ZLEVELS.size

    # 출력 Dataset (coords 먼저 명시)
    ds_out = xr.Dataset(
        coords={
            time_name: (time_name, times),
            "z": ("z", ZLEVELS, {"units": "m", "positive": "up"})
        }
    )

    # ---- RHO VARS (temp, salt)
    for vn in VARS_RHO:
        if vn not in ds:
            raise KeyError(f"'{vn}' 없음")
        v = ds[vn]
        if v.ndim != 4:
            raise ValueError(f"{vn} shape 기대: (time,s_rho,eta_rho,xi_rho), 실제: {v.shape}")

        eta_name = v.dims[-2]
        xi_name  = v.dims[-1]

        out = np.empty((Nt, Nz, v.shape[-2], v.shape[-1]), dtype=float)

        for it, t in enumerate(t_indices):
            zeta2d = zeta_all[it, :, :]
            z_rho = _make_z(ds, h, zeta2d, N, igrid=1)       # rho grid (N,eta_rho,xi_rho)
            v3d = v.isel({time_name: (t if TIME_INDEX is None else int(TIME_INDEX))}).values
            # v3d: (s_rho,eta,xi)
            out[it] = _interp_one_time(z_rho, v3d)

        ds_out[vn] = xr.DataArray(
            out,
            dims=(time_name, "z", eta_name, xi_name),
            attrs=dict(v.attrs)
        )

    # ---- U VAR
    if VAR_U in ds:
        u = ds[VAR_U]
        if u.ndim != 4:
            raise ValueError(f"{VAR_U} shape 기대: (time,s_rho,eta_u,xi_u), 실제: {u.shape}")

        eta_u = u.dims[-2]
        xi_u  = u.dims[-1]

        out_u = np.empty((Nt, Nz, u.shape[-2], u.shape[-1]), dtype=float)

        for it, t in enumerate(t_indices):
            zeta2d = zeta_all[it, :, :]
            z_u = _make_z(ds, h, zeta2d, N, igrid=3)  # u grid
            u3d = u.isel({time_name: (t if TIME_INDEX is None else int(TIME_INDEX))}).values
            if z_u.shape != u3d.shape:
                raise ValueError(f"shape mismatch u: z_u{z_u.shape} vs u3d{u3d.shape}")
            out_u[it] = _interp_one_time(z_u, u3d)

        ds_out[VAR_U] = xr.DataArray(out_u, dims=(time_name, "z", eta_u, xi_u), attrs=dict(u.attrs))

    # ---- V VAR
    if VAR_V in ds:
        v = ds[VAR_V]
        if v.ndim != 4:
            raise ValueError(f"{VAR_V} shape 기대: (time,s_rho,eta_v,xi_v), 실제: {v.shape}")

        eta_v = v.dims[-2]
        xi_v  = v.dims[-1]

        out_v = np.empty((Nt, Nz, v.shape[-2], v.shape[-1]), dtype=float)

        for it, t in enumerate(t_indices):
            zeta2d = zeta_all[it, :, :]
            z_v = _make_z(ds, h, zeta2d, N, igrid=4)  # v grid
            v3d = v.isel({time_name: (t if TIME_INDEX is None else int(TIME_INDEX))}).values
            if z_v.shape != v3d.shape:
                raise ValueError(f"shape mismatch v: z_v{z_v.shape} vs v3d{v3d.shape}")
            out_v[it] = _interp_one_time(z_v, v3d)

        ds_out[VAR_V] = xr.DataArray(out_v, dims=(time_name, "z", eta_v, xi_v), attrs=dict(v.attrs))

    # ---- 그리드/좌표 복사 (MergeError 방지: sanitize해서 넣음)
    for name in ["lon_rho", "lat_rho", "mask_rho", "lon_u", "lat_u", "lon_v", "lat_v", "mask_u", "mask_v"]:
        if name in ds:
            ds_out[name] = _coord_sanitized(ds[name])

    # ---- 참고로 h, zeta도 저장 (검증용)
    ds_out["h"] = _coord_sanitized(ds["h"])
    if TIME_INDEX is None:
        ds_out["zeta"] = xr.DataArray(zeta_all, dims=(time_name, ds["zeta"].dims[-2], ds["zeta"].dims[-1]),
                                      attrs=dict(ds["zeta"].attrs))
    else:
        ds_out["zeta"] = xr.DataArray(zeta_all[0], dims=(ds["zeta"].dims[-2], ds["zeta"].dims[-1]),
                                      attrs=dict(ds["zeta"].attrs))

    ds_out.to_netcdf(OUTPUT_NC)
    ds.close()
    print(f"written: {OUTPUT_NC}")


if __name__ == "__main__":
    main()
