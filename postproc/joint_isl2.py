# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append("C:/Users/ust21/shjo/romsforge/libs/")

import matplotlib.pyplot as plt

import utils as ut
import create_I as cn
import post_utils as pu
from scipy.interpolate import LinearNDInterpolator, griddata

pth="D:/shjo/ISL_TEST/"
pthG_src=pth+"NWP4_grd_3_10m_LP.nc"
pthG_dst=pth+"NWP12_grd_NWP4.nc"
pthI_src=pth+"NWP4_ini_DA_8401.nc"
# pthI_dst=pth+"NWP4_ini_DA_8401.nc"

pth="D:/shjo/ISL_TEST/"
pthG_src=pth+"NWP12_grd_NWP4.nc"
pthG_dst=pth+"NWP4_grd_3_10m_LP.nc"

pthI_src=pth+"NWP12_rst_8401.nc"
# pthI_src=pth+"NWP4_ini_DA_8401.nc"


cfg  = ut.parse_config("./config_all.yaml")
grd  = ut.load_roms_grid(pthG_dst)

# Read grid
ncG_src=Dataset(pthG_src)
ncG_dst=Dataset(pthG_dst)
ncI_src=Dataset(pthI_src)
# ncI_dst=Dataset(pthI_dst)

lon_src,lat_src=ncG_src['lon_rho'][:], ncG_src['lat_rho'][:]
lon_dst,lat_dst=ncG_dst['lon_rho'][:], ncG_dst['lat_rho'][:]

topo_src = ncG_src["h"][:]
topo_dst = ncG_dst["h"][:]
Rmask = ncG_dst["mask_rho"][:]
Dmask = ncG_src["mask_rho"][:]

angle_src = ncG_src["angle"][:]
angle_dst = ncG_dst["angle"][:]

vtrs=float(ncI_src["Vtransform"][:].data)
vstr=float(ncI_src["Vstretching"][:].data)
theta_s=float(ncI_src["theta_s"][:].data)
theta_b=float(ncI_src["theta_b"][:].data)
theta_b=0.1
tcline=float(ncI_src["Tcline"][:].data)
# nlayer=len(ncI_src["sc_r"][:])
nlayer=len(ncI_src["s_rho"][:])


zr_src=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer,1,topo_src,np.zeros_like(topo_src))
zr_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,20,1,topo_dst,np.zeros_like(topo_dst))
zw_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,20,5,topo_dst,np.zeros_like(topo_dst))

# --- Load fields (example: single time) ---
zeta_diff=ncI_src["zeta"][0]
temp_diff=ncI_src["temp"][0]
salt_diff=ncI_src["salt"][0]
u_diff=ncI_src["u"][0]
v_diff=ncI_src["v"][0]
ubar_diff=ncI_src["ubar"][0]
vbar_diff=ncI_src["vbar"][0]

u_diff=pu.uv2rho_rutgers_safenan(u_diff,"u")
v_diff=pu.uv2rho_rutgers_safenan(v_diff,"v")
ubar_diff=pu.uv2rho_rutgers_safenan(ubar_diff,"u")
vbar_diff=pu.uv2rho_rutgers_safenan(vbar_diff,"v")

# --- Load & calc increments ---
# factor=1
# zeta_diff=(ncI_src["zeta"][0]-ncI_src["zeta"][1])*factor
# temp_diff=(ncI_src["temp"][0]-ncI_src["temp"][1])*factor
# salt_diff=(ncI_src["salt"][0]-ncI_src["salt"][1])*factor
# u_diff=(ncI_src["u"][0]-ncI_src["u"][1])*factor
# v_diff=(ncI_src["v"][0]-ncI_src["v"][1])*factor
# ubar_diff=(ncI_src["ubar"][0]-ncI_src["ubar"][1])*factor
# vbar_diff=(ncI_src["vbar"][0]-ncI_src["vbar"][1])*factor

# u_diff=pu.uv2rho_rutgers_safenan(u_diff,"u")
# v_diff=pu.uv2rho_rutgers_safenan(v_diff,"v")
# ubar_diff=pu.uv2rho_rutgers_safenan(ubar_diff,"u")
# vbar_diff=pu.uv2rho_rutgers_safenan(vbar_diff,"v")


def _interp2d_linear_then_nearest(VD, XD, YD, XR, YR,
                                  Dmask=None,
                                  treat_zero_as_fill=False,
                                  fill_value_from_attr=None):
    x = XD.ravel(); y = YD.ravel()
    v = VD.ravel()
    if Dmask is not None:
        m = Dmask.ravel()
        valid = (m >= 0.5)
    else:
        valid = np.ones_like(v, dtype=bool)
    valid &= ~np.isnan(v)
    if fill_value_from_attr is not None:
        valid &= (v != fill_value_from_attr)
    if treat_zero_as_fill:
        valid &= (v != 0.0)
    if valid.sum() == 0:
        return np.full_like(XR, np.nan, dtype=VD.dtype)

    points = np.column_stack((x[valid], y[valid]))
    values = v[valid]
    out = griddata(points, values, (XR, YR), method='linear')
    nan_mask = np.isnan(out)
    if np.any(nan_mask):
        out[nan_mask] = griddata(points, values,
                                 (XR[nan_mask], YR[nan_mask]),
                                 method='nearest')
    return out

def interp2d_with_masks(VD, XD, YD, XR, YR, Dmask=None, Rmask=None,
                        treat_zero_as_fill=False, fill_value_from_attr=None):
    V = _interp2d_linear_then_nearest(
            VD, XD, YD, XR, YR,
            Dmask=Dmask,
            treat_zero_as_fill=treat_zero_as_fill,
            fill_value_from_attr=fill_value_from_attr)
    if Rmask is not None:
        V = np.where(Rmask >= 0.5, V, 0.0)
    return V

def interp_angle_like_matlab(angle_src, lon_src, lat_src, lon_dst, lat_dst, Dmask=None):
    valid = ~np.isnan(angle_src)
    if Dmask is not None:
        valid &= (Dmask >= 0.5)
    if valid.sum() == 0:
        return np.full_like(lon_dst, np.nan, dtype=angle_src.dtype)

    interp_lin = LinearNDInterpolator(
        (lon_src[valid], lat_src[valid]),
        angle_src[valid],
        fill_value=np.nan
    )
    ang_dst = interp_lin(lon_dst, lat_dst)
    nan_mask = np.isnan(ang_dst)
    if np.any(nan_mask):
        points  = np.column_stack((lon_src[valid], lat_src[valid]))
        values  = angle_src[valid]
        targets = (lon_dst[nan_mask], lat_dst[nan_mask])
        ang_dst[nan_mask] = griddata(points, values, targets, method='nearest')
    return ang_dst.astype(angle_src.dtype, copy=False)

def hz_to_u_v(Hz: np.ndarray):
    dzu = 0.5 * (Hz[:, :, :-1] + Hz[:, :, 1:])
    dzv = 0.5 * (Hz[:, :-1, :] + Hz[:, 1:, :])
    return dzu, dzv

def uv_barotropic_from_3d(u: np.ndarray,
                          v: np.ndarray,
                          Hz: np.ndarray,
                          mask_u: np.ndarray | None = None,
                          mask_v: np.ndarray | None = None,
                          eps: float = 1e-12):
    dzu, dzv = hz_to_u_v(Hz)
    Du = np.sum(dzu, axis=0)
    Dv = np.sum(dzv, axis=0)
    num_u = np.sum(u * dzu, axis=0)
    num_v = np.sum(v * dzv, axis=0)
    ubar = np.where(Du > eps, num_u / Du, 0.0)
    vbar = np.where(Dv > eps, num_v / Dv, 0.0)
    if mask_u is not None:
        ubar = np.where(mask_u == 1, ubar, 0.0)
    if mask_v is not None:
        vbar = np.where(mask_v == 1, vbar, 0.0)
    return ubar, vbar

# --- define field ---
field = ut.ConfigObject(zeta=zeta_diff, temp=temp_diff, salt=salt_diff,
                        u=u_diff, v=v_diff, ubar=None, vbar=None)

# --- Horizontal interp (2D: Rmask 적용, 3D: Rmask 적용하지 않음) ---
for name in ["zeta","temp","salt","u","v"]:
    if not hasattr(field, name): 
        continue
    arr = getattr(field, name)
    if arr is None:
        continue

    if arr.ndim == 2:
        out = interp2d_with_masks(arr, lon_src, lat_src, lon_dst, lat_dst,
                                  Dmask=Dmask, Rmask=Rmask,
                                  treat_zero_as_fill=True)
        setattr(field, name, out.astype(arr.dtype, copy=False))

    elif arr.ndim == 3:
        nz = arr.shape[0]
        out = np.empty((nz,) + lon_dst.shape, dtype=arr.dtype)
        for k in range(nz):
            out[k] = interp2d_with_masks(arr[k], lon_src, lat_src, lon_dst, lat_dst,
                                         Dmask=Dmask, Rmask=None,          # <-- no Rmask here
                                         treat_zero_as_fill=True)
        setattr(field, name, out)

# --- zr (donor z_r) -> dst 수평보간 (Rmask 적용하지 않음) ---
nz = zr_src.shape[0]
zr_src_remapped = np.empty((nz,) + lon_dst.shape, dtype=zr_src.dtype)
for k in range(nz):
    zr_src_remapped[k] = interp2d_with_masks(
        zr_src[k], lon_src, lat_src, lon_dst, lat_dst,
        Dmask=Dmask, Rmask=None, treat_zero_as_fill=False
    )

# --- Vertical interpolation to zr_dst ---
field2 = ut.ConfigObject()
for name in ["temp","salt","u","v"]:
    if hasattr(field, name):
        var = getattr(field, name)
        if var is None or var.ndim != 3:
            continue
        var_z = pu.vertical_interp_to_ZR(zr_src_remapped, var, zr_dst,
                                         n_jobs=-1, dedup="mean", extrap_mode="leading")
        setattr(field2, name, var_z)

# 수직보간 끝난 뒤 Rmask 적용 (MATLAB 3D 분기와 동일)
for name in ["temp","salt","u","v"]:
    if hasattr(field2, name):
        A = getattr(field2, name)
        A[:, Rmask == 0] = 0.0
        setattr(field2, name, A)

# 2D zeta는 앞 단계 결과 사용
field2.zeta = field.zeta

# --- angle: parent_angle on dst (Dmask 사용) ---
angle_par_on_dst = interp_angle_like_matlab(angle_src, lon_src, lat_src, lon_dst, lat_dst, Dmask=Dmask)

# --- Rotate vectors at rho, then stagger to U/V, then mask_u/v ---
u, v = field2.u, field2.v
Urho, Vrho = ut.rotate_vector_euler(u, v, angle_par_on_dst, to_geo=True)
Ugeo, Vgeo = ut.rotate_vector_euler(Urho, Vrho, angle_dst, to_geo=False)

mask_u = ncG_dst["mask_u"][:]
mask_v = ncG_dst["mask_v"][:]
u_on_u = ut.rho2uv(Ugeo, pos="u")
v_on_v = ut.rho2uv(Vgeo, pos="v")
u_on_u = np.where(mask_u == 1, u_on_u, 0.0)
v_on_v = np.where(mask_v == 1, v_on_v, 0.0)

field2.u = u_on_u
field2.v = v_on_v

# --- Barotropic from 3D u,v and Hz ---
Hz_dst = np.diff(zw_dst, axis=0)
ubar, vbar = uv_barotropic_from_3d(field2.u, field2.v, Hz_dst, mask_u=mask_u, mask_v=mask_v)
field2.ubar = ubar
field2.vbar  = vbar

# --- Create initial NetCDF and write ---
status = cn.create_ini__(cfg, grd, 0, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
if status:
    raise RuntimeError(f"Failed creating file {cfg.ininame}")

with Dataset(cfg.ininame, mode='a') as nc:
    nc['zeta'][0] = field2.zeta
    nc['temp'][0] = field2.temp
    nc['salt'][0] = field2.salt
    nc['u'][0]    = field2.u
    nc['v'][0]    = field2.v
    nc['ubar'][0] = field2.ubar
    nc['vbar'][0] = field2.vbar
