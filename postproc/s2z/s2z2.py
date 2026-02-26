
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ROMS s-level -> fixed z-level extraction (NO padding fill)
- boundary/out-of-range stays NaN (no surface/bottom filling)
- robust to zlevs output axis order
- robust u/v grid sizing (eta_u/xi_u, eta_v/xi_v)
"""

import os
import sys
import argparse
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
from joblib import Parallel, delayed

# romsforge libs
sys.path.append("/home/shjo/github/romsforge/libs")
import utils as ut  # zlevs()


# ============================================================
# helpers
# ============================================================
def _ensure_z_first(depth3d: np.ndarray, Ns: int) -> np.ndarray:
    """
    depth3d가 (Ns, Ny, Nx) 또는 (Ny, Nx, Ns) 등으로 올 수 있어서
    Ns가 어느 축인지 보고 (Ns, Ny, Nx)로 맞춰줌.
    """
    if depth3d.ndim != 3:
        raise ValueError(f"depth3d must be 3D, got {depth3d.ndim}D")

    if depth3d.shape[0] == Ns:
        return depth3d
    if depth3d.shape[-1] == Ns:
        return np.moveaxis(depth3d, -1, 0)
    if depth3d.shape[1] == Ns:
        # (Ny, Ns, Nx) 같은 이상 케이스
        return np.moveaxis(depth3d, 1, 0)

    raise ValueError(f"cannot find vertical axis Ns={Ns} in depth3d shape={depth3d.shape}")


def _enforce_monotonic(z: np.ndarray, v: np.ndarray, dedup="jitter", eps=1e-10):
    """
    NaN 제거 -> z 오름차순 정렬 -> 중복 처리
    (post_utils의 로직을 최소 버전으로 내장)
    """
    m = np.isfinite(z) & np.isfinite(v)
    if m.sum() < 2:
        return None, None
    z = z[m].astype(float)
    v = v[m].astype(float)

    order = np.argsort(z)
    z = z[order]
    v = v[order]

    if dedup == "jitter":
        dz = np.diff(z)
        bad = np.where(dz <= 0)[0]
        if bad.size:
            scale = max(abs(z[-1] - z[0]), 1.0)
            step = eps * scale
            for i in bad:
                if z[i + 1] <= z[i]:
                    z[i + 1] = z[i] + step
    elif dedup == "mean":
        zu, inv = np.unique(z, return_inverse=True)
        if len(zu) != len(z):
            acc = np.zeros_like(zu, dtype=float)
            cnt = np.zeros_like(zu, dtype=int)
            for t, val in zip(inv, v):
                acc[t] += val
                cnt[t] += 1
            z = zu
            v = acc / np.maximum(cnt, 1)

    if z.size < 2:
        return None, None
    return z, v


def interpolate_s_to_zlevels_nan(
    depth3d: np.ndarray,   # (Ns, Ny, Nx)
    roms_var: np.ndarray,  # (Ns, Ny, Nx)
    stdvdepth: np.ndarray, # (Nz,)
    *,
    n_jobs: int = -1,
    dedup: str = "jitter"
) -> np.ndarray:
    """
    s->z 보간. out-of-range는 끝까지 NaN (패딩/표층저층 채움 없음).
    """
    depth3d = np.asarray(depth3d)
    roms_var = np.asarray(roms_var)
    stdvdepth = np.asarray(stdvdepth, dtype=float)

    if depth3d.shape != roms_var.shape:
        raise ValueError(f"depth3d shape {depth3d.shape} != roms_var shape {roms_var.shape}")

    Ns, Ny, Nx = depth3d.shape
    Nz = stdvdepth.size

    # stdvdepth는 어떤 순서든 가능하게: 내부는 오름차순으로 계산 후 복원
    if Nz == 0:
        return np.empty((0, Ny, Nx), dtype=float)

    if np.all(np.diff(stdvdepth) > 0):
        std_sorted = stdvdepth
        restore = None
    else:
        order = np.argsort(stdvdepth)
        std_sorted = stdvdepth[order]
        restore = np.argsort(order)

    def interp_one_x(ix: int) -> np.ndarray:
        out = np.full((Nz, Ny), np.nan, dtype=float)
        for jy in range(Ny):
            z_col = depth3d[:, jy, ix]
            v_col = roms_var[:, jy, ix]

            z, v = _enforce_monotonic(z_col, v_col, dedup=dedup)
            if z is None:
                continue

            # out-of-range는 NaN 유지
            f = interp1d(
                z, v, kind="linear",
                bounds_error=False,
                fill_value=(np.nan, np.nan),
                assume_sorted=True
            )
            col = f(std_sorted)

            # 원 stdvdepth 순서 복원
            if restore is not None:
                tmp = np.empty_like(col)
                tmp[restore] = col
                col = tmp

            out[:, jy] = col
        return out

    cols = Parallel(n_jobs=n_jobs, prefer="threads")(delayed(interp_one_x)(i) for i in range(Nx))
    out = np.stack(cols, axis=-1)  # (Nz, Ny, Nx)
    return out


def _avg_xi(a):
    # rho -> u (xi-1)
    return 0.5 * (a[..., :, 1:] + a[..., :, :-1])


def _avg_eta(a):
    # rho -> v (eta-1)
    return 0.5 * (a[..., 1:, :] + a[..., :-1, :])


def _get_time_name(ds: xr.Dataset) -> str:
    for cand in ["ocean_time", "time"]:
        if cand in ds.dims:
            return cand
        if cand in ds.coords:
            return cand
    raise KeyError("time dim/coord not found (ocean_time/time)")


def _pick_zlevels(zlevels_arg: str) -> np.ndarray:
    """
    --zlevels "-5,-10,-20" or "--zlevels_file z.txt" 같은 식으로 확장 가능하게
    일단은 콤마 리스트만 받자.
    """
    parts = [p.strip() for p in zlevels_arg.split(",") if p.strip()]
    z = np.array([float(p) for p in parts], dtype=float)
    if z.size == 0:
        raise ValueError("zlevels empty")
    return z


# ============================================================
# main
# ============================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--his", required=True, help="ROMS history file path")
    ap.add_argument("--out", required=True, help="output NetCDF path")
    ap.add_argument("--zlevels", required=True, help='comma list, e.g. "-5,-10,-20,-50" (meters, negative down)')
    ap.add_argument("--vars", default="temp,salt,u,v,zeta", help="comma list of variables to export")
    ap.add_argument("--tidx", default=None, help="time indices like '0' or '0:10' or '0,3,5' (default: all)")
    ap.add_argument("--n_jobs", type=int, default=-1)
    args = ap.parse_args()

    his = args.his
    out_nc = args.out
    zlevels = _pick_zlevels(args.zlevels)
    varlist = [v.strip() for v in args.vars.split(",") if v.strip()]

    ds = xr.open_dataset(his, decode_times=False)
    tname = _get_time_name(ds)

    # dims
    if "s_rho" not in ds.dims:
        raise KeyError("dataset has no s_rho dim")
    Ns = int(ds.sizes["s_rho"])

    # required grid vars
    for req in ["h", "zeta", "Vtransform", "Vstretching", "theta_s", "theta_b", "hc"]:
        if req not in ds:
            raise KeyError(f"missing required ROMS var: {req}")

    # base grid
    h_rho = ds["h"].values
    eta_rho, xi_rho = h_rho.shape

    # zeta time-dependent
    # ROMS params
    Vtransform = int(np.array(ds["Vtransform"].values).ravel()[0])
    Vstretching = int(np.array(ds["Vstretching"].values).ravel()[0])
    theta_s = float(np.array(ds["theta_s"].values).ravel()[0])
    theta_b = float(np.array(ds["theta_b"].values).ravel()[0])
    hc = float(np.array(ds["hc"].values).ravel()[0])

    # coords (keep as normal data vars; don't rely on xarray coords merge)
    lon_rho = ds["lon_rho"].values if "lon_rho" in ds else None
    lat_rho = ds["lat_rho"].values if "lat_rho" in ds else None
    lon_u = ds["lon_u"].values if "lon_u" in ds else None
    lat_u = ds["lat_u"].values if "lat_u" in ds else None
    lon_v = ds["lon_v"].values if "lon_v" in ds else None
    lat_v = ds["lat_v"].values if "lat_v" in ds else None

    # time indices selection
    Nt = int(ds.sizes[tname])
    if args.tidx is None:
        tidxs = list(range(Nt))
    else:
        s = args.tidx.strip()
        if ":" in s:
            a, b = s.split(":")
            a = int(a) if a else 0
            b = int(b) if b else Nt
            tidxs = list(range(a, min(b, Nt)))
        elif "," in s:
            tidxs = [int(x) for x in s.split(",") if x.strip() != ""]
            tidxs = [x for x in tidxs if 0 <= x < Nt]
        else:
            tidxs = [int(s)]
    if len(tidxs) == 0:
        raise ValueError("no valid time indices selected")

    # allocate output arrays
    Nz = len(zlevels)
    out_vars = {}

    # precompute static h_u, h_v
    if "u" in ds:
        eta_u = int(ds["u"].sizes.get("eta_u", ds["u"].shape[-2]))
        xi_u = int(ds["u"].sizes.get("xi_u", ds["u"].shape[-1]))
        h_u = _avg_xi(h_rho)[:eta_u, :xi_u]
    else:
        h_u = None

    if "v" in ds:
        eta_v = int(ds["v"].sizes.get("eta_v", ds["v"].shape[-2]))
        xi_v = int(ds["v"].sizes.get("xi_v", ds["v"].shape[-1]))
        h_v = _avg_eta(h_rho)[:eta_v, :xi_v]
    else:
        h_v = None

    # output dims names (preserve if present)
    eta_rho_name = "eta_rho" if "eta_rho" in ds.dims else "eta"
    xi_rho_name  = "xi_rho"  if "xi_rho"  in ds.dims else "xi"

    # Build output dataset explicitly (avoid merge hell)
    ds_out = xr.Dataset()
    ds_out = ds_out.assign_coords({tname: ds[tname].values, "z": zlevels})

    # attach lon/lat as data vars (not coords) to avoid coordinate merge ambiguity
    if lon_rho is not None and lat_rho is not None:
        ds_out["lon_rho"] = xr.DataArray(lon_rho, dims=(eta_rho_name, xi_rho_name))
        ds_out["lat_rho"] = xr.DataArray(lat_rho, dims=(eta_rho_name, xi_rho_name))
    if lon_u is not None and lat_u is not None and ("u" in ds):
        ds_out["lon_u"] = xr.DataArray(lon_u, dims=("eta_u", "xi_u"))
        ds_out["lat_u"] = xr.DataArray(lat_u, dims=("eta_u", "xi_u"))
    if lon_v is not None and lat_v is not None and ("v" in ds):
        ds_out["lon_v"] = xr.DataArray(lon_v, dims=("eta_v", "xi_v"))
        ds_out["lat_v"] = xr.DataArray(lat_v, dims=("eta_v", "xi_v"))

    # copy a few scalars safely
    for name in ["Vtransform", "Vstretching", "theta_s", "theta_b", "hc"]:
        if name in ds:
            ds_out[name] = ds[name]

    # zeta (optional)
    if "zeta" in varlist and "zeta" in ds:
        ds_out["zeta"] = xr.DataArray(ds["zeta"].values, dims=(tname, eta_rho_name, xi_rho_name), attrs=ds["zeta"].attrs)

    # process each time index
    # We'll store outputs in full time length, but only fill selected tidxs
    def _alloc_4d(dims):
        # dims = (Nt, Nz, Ny, Nx)
        arr = np.full(dims, np.nan, dtype=np.float32)
        return arr

    # rho-grid 3D vars
    rho3d_names = [vn for vn in varlist if vn in ds and vn not in ("u", "v", "zeta")]
    # but keep only those that have s_rho dim
    rho3d_names = [vn for vn in rho3d_names if ("s_rho" in ds[vn].dims)]

    # allocate
    for vn in rho3d_names:
        Ny = int(ds[vn].sizes.get(eta_rho_name, ds[vn].shape[-2]))
        Nx = int(ds[vn].sizes.get(xi_rho_name, ds[vn].shape[-1]))
        out_vars[vn] = _alloc_4d((Nt, Nz, Ny, Nx))

    if "u" in varlist and "u" in ds:
        out_u = _alloc_4d((Nt, Nz, eta_u, xi_u))
    else:
        out_u = None

    if "v" in varlist and "v" in ds:
        out_v = _alloc_4d((Nt, Nz, eta_v, xi_v))
    else:
        out_v = None

    # interpolate per-time
    for it in tidxs:
        # zeta_rho at time it
        zeta_rho = ds["zeta"].isel({tname: it}).values

        # --- rho depth ---
        z_rho = ut.zlevs(Vtransform, Vstretching, theta_s, theta_b, hc, Ns, 1, h_rho, zeta_rho)
        z_rho = _ensure_z_first(np.asarray(z_rho), Ns)  # (Ns, eta_rho, xi_rho)

        # rho vars
        for vn in rho3d_names:
            v3d = ds[vn].isel({tname: it}).values  # (Ns, eta, xi) in most files
            v3d = _ensure_z_first(np.asarray(v3d), Ns)
            out = interpolate_s_to_zlevels_nan(z_rho, v3d, zlevels, n_jobs=args.n_jobs)
            out_vars[vn][it, :, :, :] = out.astype(np.float32)

        # --- u ---
        if out_u is not None:
            # zeta_u, depth_u
            zeta_u = _avg_xi(zeta_rho)[:eta_u, :xi_u]
            z_u = ut.zlevs(Vtransform, Vstretching, theta_s, theta_b, hc, Ns, 2, h_u, zeta_u)
            z_u = _ensure_z_first(np.asarray(z_u), Ns)  # (Ns, eta_u, xi_u)

            u3d = ds["u"].isel({tname: it}).values
            u3d = _ensure_z_first(np.asarray(u3d), Ns)
            # just in case: slice to match exactly
            u3d = u3d[:, :eta_u, :xi_u]
            z_u = z_u[:, :eta_u, :xi_u]

            out = interpolate_s_to_zlevels_nan(z_u, u3d, zlevels, n_jobs=args.n_jobs)
            out_u[it, :, :, :] = out.astype(np.float32)

        # --- v ---
        if out_v is not None:
            zeta_v = _avg_eta(zeta_rho)[:eta_v, :xi_v]
            z_v = ut.zlevs(Vtransform, Vstretching, theta_s, theta_b, hc, Ns, 3, h_v, zeta_v)
            z_v = _ensure_z_first(np.asarray(z_v), Ns)  # (Ns, eta_v, xi_v)

            v3d = ds["v"].isel({tname: it}).values
            v3d = _ensure_z_first(np.asarray(v3d), Ns)
            v3d = v3d[:, :eta_v, :xi_v]
            z_v = z_v[:, :eta_v, :xi_v]

            out = interpolate_s_to_zlevels_nan(z_v, v3d, zlevels, n_jobs=args.n_jobs)
            out_v[it, :, :, :] = out.astype(np.float32)

        print(f"[OK] time index {it}/{Nt-1}")

    # write outputs to ds_out
    for vn, arr in out_vars.items():
        ds_out[vn] = xr.DataArray(arr, dims=(tname, "z", eta_rho_name, xi_rho_name), attrs=ds[vn].attrs)

    if out_u is not None:
        ds_out["u"] = xr.DataArray(out_u, dims=(tname, "z", "eta_u", "xi_u"), attrs=ds["u"].attrs)

    if out_v is not None:
        ds_out["v"] = xr.DataArray(out_v, dims=(tname, "z", "eta_v", "xi_v"), attrs=ds["v"].attrs)

    # save
    os.makedirs(os.path.dirname(out_nc) or ".", exist_ok=True)
    ds_out.to_netcdf(out_nc)
    ds.close()
    ds_out.close()

    print("[saved]", out_nc)


if __name__ == "__main__":
    main()
