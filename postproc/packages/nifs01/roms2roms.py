# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))

import matplotlib.pyplot as plt

import utils as ut
import create_I as cn
import post_utils as pu
from scipy.interpolate import LinearNDInterpolator, griddata

cfg  = ut.parse_config("./config_isl_h2l.yaml")
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
# idt=0
TIME=Dataset(cfg["ininame_src"])["ocean_time"]
ocean_time=TIME[:]
time_ref=TIME.units
with Dataset(cfg["ininame_src"], maskandscale=True) as nc_raw:
    nc = ut.MaskedNetCDF(nc_raw)
    zeta = nc.get('zeta', slice(None)).squeeze()
    ubar = nc.get('ubar', slice(None)).squeeze()
    vbar = nc.get('vbar', slice(None)).squeeze()
    temp = nc.get('temp', slice(None)).squeeze()
    salt = nc.get('salt', slice(None)).squeeze()
    u    = nc.get('u',    slice(None)).squeeze()
    v    = nc.get('v',    slice(None)).squeeze()

    # U/V -> RHO (NaN 보존)
    u=pu.uv2rho_rutgers_safenan(u,"u")
    v=pu.uv2rho_rutgers_safenan(v,"v")
    ubar=pu.uv2rho_rutgers_safenan(ubar,"u")
    vbar=pu.uv2rho_rutgers_safenan(vbar,"v")

    field = ut.ConfigObject(zeta=zeta, ubar=None, vbar=None, temp=temp, salt=salt, u=u, v=v)


# --- 캐시 생성 (한 번만) ---
cache_rho = pu.build_lni_cache(lon_src, lat_src, Dmask, lon_dst, lat_dst)

# --- Horizontal interp (2D: Rmask 적용, 3D: Rmask 적용하지 않음) ---
for name in ["zeta","temp","salt","u","v"]:
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
for name in ["temp","salt","u","v"]:
    print('Vertical intrp: '+name)
    if hasattr(field, name):
        var = getattr(field, name)
        if var is None or var.ndim != 3:
            print(name)
            continue
        # MATLAB 경계복제와 유사한 padding 외삽 사용(너가 쓰던 설정 유지)
        var_z = pu.vertical_interp_to_ZR(
            zr_src_remapped.astype(np.float64, copy=False),
            var.astype(np.float64, copy=False),
            zr_dst.astype(np.float64, copy=False),
            n_jobs=-1, dedup="mean", extrap_mode="padding",
            zsur=0.0, zbot=None
        )
        # var_z = pu.vertical_interp_to_ZR(zr_src_remapped, var, zr_dst,
        #                          n_jobs=-1, dedup="mean", extrap_mode="leading")
        setattr(field2, name, var_z.astype(np.float64, copy=False))

# 수직보간 끝난 뒤 Rmask 적용 (스칼라만)
for name in ["temp","salt"]:
    if hasattr(field2, name):
        A = getattr(field2, name)
        A[:, Rmask == 0] = 0.0
        setattr(field2, name, A)

# 2D zeta는 앞 단계 결과 사용
field2.zeta = field.zeta

# --- angle: parent_angle on dst (Dmask 사용; 캐시 재사용) ---
angle_par_on_dst = pu.apply_lni_cache(angle_src.astype(np.float64, copy=False), cache_rho)

# --- Rotate vectors at rho, then stagger to U/V, then mask_u/v ---
u, v = field2.u, field2.v
Urho, Vrho = ut.rotate_vector_euler(u, v, angle_par_on_dst, to_geo=True)
Ugeo, Vgeo = ut.rotate_vector_euler(Urho, Vrho, angle_dst, to_geo=False)

mask_u = grd_dst["mask_u"][:]
mask_v = grd_dst["mask_v"][:]

u_on_u = pu.rho2uv_nanaware(Ugeo, pos="u")
v_on_v = pu.rho2uv_nanaware(Vgeo, pos="v")

# 마지막에만 mask_u/v 적용
u_on_u = np.where(mask_u == 1, u_on_u, 0.0)
v_on_v = np.where(mask_v == 1, v_on_v, 0.0)

field2.u = u_on_u.astype(np.float64, copy=False)
field2.v = v_on_v.astype(np.float64, copy=False)

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
    for varname in vars(field):
        nc[varname][0] = field2[varname]



















