# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 11:20:08 2025

@author: ust21
"""

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


cfg  = ut.parse_config("C:/Users/ust21/shjo/romsforge/postproc/ISL_system/config_ISL.yaml")
grd_src  = ut.load_roms_grid(cfg["grdname_src"])
grd_dst  = ut.load_roms_grid(cfg["grdname_dst"])
# ncI_dst=Dataset(pthI_dst)

# Read grid
ncI_src=Dataset(cfg["ininame_src"])

lon_src,lat_src=grd_src['lon'][:], grd_src['lat'][:]
lon_dst,lat_dst=grd_dst['lon'][:], grd_dst['lat'][:]

topo_src = grd_src["topo"][:]
topo_dst = grd_dst["topo"][:]
Rmask = grd_dst["mask_rho"][:]
Dmask = grd_src["mask_rho"][:]

angle_src = grd_src["angle"][:]
angle_dst = grd_dst["angle"][:]

vtrs=cfg["vertical_src"]["vtransform"]
vstr=cfg["vertical_src"]["vstretching"]
theta_s=cfg["vertical_src"]["theta_s"]
theta_b=cfg["vertical_src"]["theta_b"]
tcline=cfg["vertical_src"]["tcline"]
nlayer_src=cfg["vertical_src"]["layer_n"]
nlayer_dst=cfg["vertical_dst"]["layer_n"]

zr_src=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer_src,1,topo_src,np.zeros_like(topo_src))
zr_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer_dst,1,topo_dst,np.zeros_like(topo_dst))
zw_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer_dst,5,topo_dst,np.zeros_like(topo_dst))


# from scipy import io
# MAT=io.loadmat("D:/shjo/ISL_TEST/test/zr.mat")
# _zr=MAT["zr"].transpose([2,1,0])
# diff=_zr-zr_dst
# diff[diff==0]=np.nan
# plt.pcolor(diff[0],vmin=-0.1,vmax=0.1,cmap='bwr'); plt.colorbar()

# plt.pcolor(zr_dst[-1],vmin=-5,vmax=0)
# plt.pcolor(_zr[-1],vmin=-5,vmax=0)



# --- Load fields (example: single time) ---
# idt=0
with Dataset(cfg["ininame_src"], maskandscale=True) as nc_raw:
    nc = ut.MaskedNetCDF(nc_raw)
    zeta_post = nc.get('zeta', 0, slice(None))
    ubar_post = nc.get('ubar', 0, slice(None))
    vbar_post = nc.get('vbar', 0, slice(None))
    temp_post = nc.get('temp', 0, slice(None))
    salt_post = nc.get('salt', 0, slice(None))
    u_post    = nc.get('u',    0, slice(None))
    v_post    = nc.get('v',    0, slice(None))
    
    zeta_prior = nc.get('zeta', 1, slice(None))
    ubar_prior = nc.get('ubar', 1, slice(None))
    vbar_prior = nc.get('vbar', 1, slice(None))
    temp_prior = nc.get('temp', 1, slice(None))
    salt_prior = nc.get('salt', 1, slice(None))
    u_prior    = nc.get('u',    1, slice(None))
    v_prior    = nc.get('v',    1, slice(None))

    zeta = zeta_post - zeta_prior 
    ubar = ubar_post - ubar_prior 
    vbar = vbar_post - vbar_prior 
    temp = temp_post - temp_prior 
    salt = salt_post - salt_prior 
    u = u_post - u_prior 
    v = v_post - v_prior 

    # U/V -> RHO (NaN 보존)
    u=pu.uv2rho_rutgers_safenan(u,"u")
    v=pu.uv2rho_rutgers_safenan(v,"v")
    ubar=pu.uv2rho_rutgers_safenan(ubar,"u")
    vbar=pu.uv2rho_rutgers_safenan(vbar,"v")

    field = ut.ConfigObject(zeta=zeta, ubar=None, vbar=None, temp=temp, salt=salt, u=u, v=v)




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
    out = griddata(points, values, (XR, YR), method='linear',rescale=True)
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

# # --- define field ---
# field = ut.ConfigObject(zeta=zeta_diff, temp=temp_diff, salt=salt_diff,
#                         u=u_diff, v=v_diff, ubar=None, vbar=None)

# --- Horizontal interp (2D: Rmask 적용, 3D: Rmask 적용하지 않음) ---
# for name in ["zeta","temp","salt","u","v"]:
for name in ["temp"]:

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

from scipy import io
MAT=io.loadmat("D:/shjo/ISL_TEST/test/temp20_zr.mat")
_temp=MAT["V"].transpose([2,1,0])
_Z=MAT["Z"].transpose([2,1,0])

out.shape
diff=_temp-out
diff[diff==0]=np.nan
plt.pcolor(diff[-1],vmin=-0.1,vmax=0.1,cmap='bwr'); plt.colorbar()


nz = arr.shape[0]
out = np.empty((nz,) + lon_dst.shape, dtype=arr.dtype)
for k in range(nz):
    out[k] = griddata( (lon_src.flatten(), lat_src.flatten()), arr[k].flatten(), (lon_dst, lat_dst),method='linear',rescale=True)


# --- zr (donor z_r) -> dst 수평보간 (Rmask 적용하지 않음) ---
nz = zr_src.shape[0]
zr_src_remapped = np.empty((nz,) + lon_dst.shape, dtype=zr_src.dtype)
for k in range(nz):
    zr_src_remapped[k] = interp2d_with_masks(
        zr_src[k], lon_src, lat_src, lon_dst, lat_dst,
        Dmask=Dmask, Rmask=None, treat_zero_as_fill=False
    )
zr_src_remapped[-1]=0
    
#
# out = np.empty((nz,) + lon_dst.shape, dtype=arr.dtype)
# for k in range(nz):
#     out[k] = griddata( (lon_src.flatten(), lat_src.flatten()), zr_src[k].flatten(), (lon_dst, lat_dst))


# _Z.shape
# diff= zr_src_remapped - _Z
# plt.pcolor(diff[-10],vmin=-10,vmax=10,cmap='bwr'); plt.colorbar()

# plt.pcolor(out[-1])
# plt.pcolor(_Z[-1])

    

# --- Vertical interpolation to zr_dst ---
field2 = ut.ConfigObject()
for name in ["temp"]:
    if hasattr(field, name):
        var = getattr(field, name)
        if var is None or var.ndim != 3:
            continue
        var_z = pu.vertical_interp_to_ZR(zr_src_remapped, var, zr_dst,
                                         n_jobs=-1, dedup="mean", extrap_mode="leading")
        setattr(field2, name, var_z)

# 수직보간 끝난 뒤 Rmask 적용 (MATLAB 3D 분기와 동일)
for name in ["temp"]:
    if hasattr(field2, name):
        A = getattr(field2, name)
        A[:, Rmask == 0] = 0.0
        setattr(field2, name, A)

A.shape
MAT=io.loadmat("D:/shjo/ISL_TEST/test/temp_intrp.mat")["V"].transpose([2,1,0])
diff=var_z-MAT
plt.pcolor(diff[0],vmin=-1,vmax=1,cmap='bwr')




# 2D zeta는 앞 단계 결과 사용
field2.zeta = field.zeta

# --- angle: parent_angle on dst (Dmask 사용) ---
angle_par_on_dst = interp_angle_like_matlab(angle_src, lon_src, lat_src, lon_dst, lat_dst, Dmask=Dmask)
# from scipy import io
# MAT=io.loadmat("D:/shjo/ISL_TEST/test/angle.mat")["parent_angle"].transpose()
# plt.pcolor(angle_par_on_dst-MAT)
# plt.colorbar()
# plt.pcolor(MAT)

# --- Rotate vectors at rho, then stagger to U/V, then mask_u/v ---
u, v = field2.u, field2.v
Urho, Vrho = ut.rotate_vector_euler(u, v, angle_par_on_dst, to_geo=True)
Ugeo, Vgeo = ut.rotate_vector_euler(Urho, Vrho, angle_dst, to_geo=False)

mask_u = grd_dst["mask_u"][:]
mask_v = grd_dst["mask_v"][:]
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
with Dataset(cfg["ininame_dst"], maskandscale=True) as nc_raw:
    nc = ut.MaskedNetCDF(nc_raw)
    zeta_hres = nc.get('zeta', 0, slice(None))
    ubar_hres = nc.get('ubar', 0, slice(None))
    vbar_hres = nc.get('vbar', 0, slice(None))
    temp_hres = nc.get('temp', 0, slice(None))
    salt_hres = nc.get('salt', 0, slice(None))
    u_hres    = nc.get('u',    0, slice(None))
    v_hres    = nc.get('v',    0, slice(None))


add_zeta_incre = zeta_hres + getattr(field2,'zeta')
setattr(field2,"zeta",add_zeta_incre)
add_ubar_incre = ubar_hres + getattr(field2,'ubar')
setattr(field2,"ubar",add_ubar_incre)
add_vbar_incre = vbar_hres + getattr(field2,'vbar')
setattr(field2,"vbar",add_vbar_incre)
add_temp_incre = temp_hres + getattr(field2,'temp')
setattr(field2,"temp",add_temp_incre)
add_salt_incre = salt_hres + getattr(field2,'salt')
setattr(field2,"salt",add_salt_incre)
add_u_incre = u_hres + getattr(field2,'u')
setattr(field2,"u",add_u_incre)
add_v_incre = v_hres + getattr(field2,'v')
setattr(field2,"v",add_v_incre)


with Dataset(cfg["ininame_dst"], mode='a') as nc:
    nc['zeta'][0] = field2.zeta.astype(nc['zeta'].dtype, copy=False)
    nc['temp'][0] = field2.temp.astype(nc['temp'].dtype, copy=False)
    nc['salt'][0] = field2.salt.astype(nc['salt'].dtype, copy=False)
    nc['u'][0]    = field2.u.astype(nc['u'].dtype, copy=False)
    nc['v'][0]    = field2.v.astype(nc['v'].dtype, copy=False)
    nc['ubar'][0] = field2.ubar.astype(nc['ubar'].dtype, copy=False)
    nc['vbar'][0] = field2.vbar.astype(nc['vbar'].dtype, copy=False)

