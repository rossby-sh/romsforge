# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
from pathlib import Path
base = Path(__file__).resolve().parent
import matplotlib.pyplot as plt

import utils as ut
import create_I as cn
import post_utils as pu
from scipy.interpolate import LinearNDInterpolator, griddata

cfg  = ut.parse_config(str(base / "config_sigma2sigma.yaml"))
grd_src  = ut.load_roms_grid(cfg["grdname_src"])
grd_dst  = ut.load_roms_grid(cfg["grdname_dst"])

# Read grid
ncI_src=Dataset(cfg["ininame_src"])

lon_src,lat_src=grd_src['lon'][:], grd_src['lat'][:]
lon_dst,lat_dst=grd_dst['lon'][:], grd_dst['lat'][:]

topo_src = grd_src["topo"][:]
topo_dst = grd_dst["topo"][:]
# topo_dst = topo_src_h # 15km -> 5km 인터폴레
Rmask = grd_dst["mask_rho"][:]
Dmask = grd_src["mask_rho"][:]

lon_dst_u = grd_dst["lon_u"][:]
lat_dst_u = grd_dst["lat_u"][:]
lon_dst_v = grd_dst["lon_v"][:]
lat_dst_v = grd_dst["lat_v"][:]

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


# zr_src=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer_src,1,topo_src,Dataset(cfg["ininame_src"])['zeta'][0])
# zr_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer_dst,1,topo_dst,field.zeta)
# zw_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer_dst,5,topo_dst,field.zeta)

# --- Load fields (example: single time) ---
idt=0
TIME=Dataset(cfg["ininame_src"])["ocean_time"]

ocean_time=TIME[idt]

time_ref=TIME.units
with Dataset(cfg["ininame_src"], maskandscale=True) as nc_raw:
    nc = ut.MaskedNetCDF(nc_raw)
    zeta = nc.get('zeta',idt, slice(None)).squeeze()
    # NOTE: ubar/vbar in source may be grid-aligned (u/v points); we ignore and recompute from 3D for consistency.
    # NOTE: ubar/vbar in source may be grid-aligned (u/v points); we ignore and recompute from 3D for consistency.
    temp = nc.get('temp',idt, slice(None)).squeeze()
    salt = nc.get('salt',idt, slice(None)).squeeze()
    u_east = nc.get('u_eastward',idt, slice(None)).squeeze()
    v_north = nc.get('v_northward',idt, slice(None)).squeeze()

    NO3 = nc.get('NO3',idt, slice(None)).squeeze()
    phyt = nc.get('phytoplankton',idt, slice(None)).squeeze()
    zoop = nc.get('zooplankton',idt, slice(None)).squeeze()
    detr = nc.get('detritus',idt, slice(None)).squeeze()

#    phytFe = nc.get('phytoplanktonFe',0, slice(None)).squeeze()
#    iron = nc.get('iron',0, slice(None)).squeeze()

    # U/V -> RHO (NaN 보존)
    # u_eastward/v_northward are already at rho-points in geographic (E/N) components; no uv2rho needed.
    # u_eastward/v_northward are already at rho-points in geographic (E/N) components; no uv2rho needed.
    # (ubar/vbar ignored; recomputed later from 3D)
    # (ubar/vbar ignored; recomputed later from 3D)

    field = ut.ConfigObject(zeta=zeta, ubar=None, vbar=None, temp=temp, salt=salt, u_eastward=u_east, v_northward=v_north,\
            phytoplankton=phyt,zooplankton=zoop,NO3=NO3,detritus=detr)


# --- 캐시 생성 (한 번만) ---
cache_rho = pu.build_lni_cache(lon_src, lat_src, Dmask, lon_dst, lat_dst)

# --- Horizontal interp (2D: Rmask 적용, 3D: Rmask 적용하지 않음) ---
for name in ["zeta","temp","salt","u_eastward","v_northward","zooplankton","phytoplankton","NO3",\
        "detritus"]:
    print("Horizontal intrp: "+name)
    if not hasattr(field, name):
        continue
    arr = getattr(field, name)
    if arr is None:
        continue

    if arr.ndim == 2:
        out = pu.apply_lni_cache(arr.astype(np.float64, copy=False), cache_rho)
        # 2D는 여기서 Rmask 적용 가능
        out = np.where(Rmask==1, out, 0.0)
        setattr(field, name, out.astype(np.float64, copy=False))

    elif arr.ndim == 3:
        nz = arr.shape[0]
        out = np.empty((nz,) + lon_dst.shape, dtype=np.float64)
        for k in range(nz):
            out[k] = pu.apply_lni_cache(arr[k].astype(np.float64, copy=False), cache_rho)
        # 3D는 아직 Rmask 적용하지 않음 (수직보간 뒤 한 번에)
        setattr(field, name, out)

# --- zr (donor z_r) -> dst 수평보간 (Rmask 적용하지 않음) ---
nz = zr_src.shape[0]
zr_src_remapped = np.empty((nz,) + lon_dst.shape, dtype=np.float64)
for k in range(nz):
    zr_src_remapped[k] = pu.apply_lni_cache(zr_src[k].astype(np.float64, copy=False), cache_rho)
# zr_src_remapped[-1]=0
# --- Vertical interpolation to zr_dst ---
field2 = ut.ConfigObject()
for name in ["temp","salt","u_eastward","v_northward","zooplankton","phytoplankton","NO3",\
        "detritus"]:
    print('Vertical intrp: '+name)
    if hasattr(field, name):
        var = getattr(field, name)
        if var is None or var.ndim != 3:
            print(var.shape)
            print(name)
            continue
        # MATLAB 경계복제와 유사한 padding 외삽 사용(너가 쓰던 설정 유지)
        #var_z = pu.vertical_interp_to_ZR(
        #    zr_src_remapped.astype(np.float64, copy=False),
        #    var.astype(np.float64, copy=False),
        #    zr_dst.astype(np.float64, copy=False),
        #    n_jobs=-1, dedup="mean", extrap_mode="padding",
        #    zsur=0.0, zbot=6000
        #)
        var_z = pu.vertical_interp_to_ZR(zr_src_remapped, var, zr_dst,n_jobs=-1, dedup="mean", extrap_mode="leading")
        setattr(field2, name, var_z.astype(np.float64, copy=False))

# 수직보간 끝난 뒤 Rmask 적용 (스칼라만)
for name in ["temp","salt","zooplankton","phytoplankton","NO3",\
        "detritus"]:
    if hasattr(field2, name):
        #print(name)
        #if name=="iron":
        #    print(getattr(field2,name))
        A = getattr(field2, name)
        A[:, Rmask == 0] = 0.0
        setattr(field2, name, A)

# 2D zeta는 앞 단계 결과 사용
field2.zeta = field.zeta

#print(vars(field))
#print(vars(field2))

# --- Convert geographic (E/N) at rho to ROMS grid-aligned (xi/eta), then stagger to U/V ---
# NOTE: ROMS expects u,v aligned with grid (xi/eta) on u/v points, not eastward/northward on rho points.



# Geographic components on rho after vertical interp
u_east_rho = field2.u_eastward
v_north_rho = field2.v_northward

# Convert (E,N) -> (xi,eta) at rho using dst grid angle (angle between xi-axis and East, CCW)
ang = angle_dst.astype(np.float64, copy=False)
u_xi_rho  =  u_east_rho*np.cos(ang) + v_north_rho*np.sin(ang)
v_eta_rho = -u_east_rho*np.sin(ang) + v_north_rho*np.cos(ang)

# Stagger to u/v points (NaN-aware) and apply masks at the end
mask_u = grd_dst['mask_u'][:]
mask_v = grd_dst['mask_v'][:]
u_on_u = pu.rho2uv_nanaware(u_xi_rho,  pos='u')
v_on_v = pu.rho2uv_nanaware(v_eta_rho, pos='v')

# 마지막에만 mask_u/v 적용
u_on_u = np.where(mask_u == 1, u_on_u, 0.0)
v_on_v = np.where(mask_v == 1, v_on_v, 0.0)

u_on_u = ut.flood_horizontal(u_on_u, lon_dst_u, lat_dst_u, 'edt')
v_on_v = ut.flood_horizontal(v_on_v, lon_dst_v, lat_dst_v, 'edt')

field2.u = u_on_u.astype(np.float64, copy=False)
field2.v = v_on_v.astype(np.float64, copy=False)



# 마지막에만 mask_u/v 적용



# --- Barotropic from 3D u,v and Hz ---
Hz_dst = np.diff(zw_dst, axis=0).astype(np.float64, copy=False)
ubar, vbar = pu.uv_barotropic_from_3d(field2.u, field2.v, Hz_dst, mask_u=mask_u, mask_v=mask_v)
field2.ubar = ubar.astype(np.float64, copy=False)
field2.vbar = vbar.astype(np.float64, copy=False)

# load hres initial
status = cn.create_roms(cfg, grd_dst, ocean_time,time_ref, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
if status:
    raise RuntimeError(f"Failed creating file {cfg.ininame}")

with Dataset(cfg["ininame_dst"], mode='a') as nc:
    for varname in vars(field2):
        if varname not in nc.variables:
            continue
        nc[varname][0] = field2[varname]



















