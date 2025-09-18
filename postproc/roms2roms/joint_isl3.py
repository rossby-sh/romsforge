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
from scipy.spatial import Delaunay, cKDTree  # ← 추가

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
nlayer=len(ncI_src["s_rho"][:])

zr_src=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer,1,topo_src,np.zeros_like(topo_src))
zr_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,20,1,topo_dst,np.zeros_like(topo_dst))
zw_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,20,5,topo_dst,np.zeros_like(topo_dst))

# --- Load fields (example: single time) ---
idt=0
with Dataset(pthI_src, maskandscale=True) as nc_raw:
    nc = ut.MaskedNetCDF(nc_raw)
    zeta = nc.get('zeta', idt, slice(None))
    ubar = nc.get('ubar', idt, slice(None))
    vbar = nc.get('vbar', idt, slice(None))

    temp = nc.get('temp', idt, slice(None))
    salt = nc.get('salt', idt, slice(None))
    u    = nc.get('u',    idt, slice(None))
    v    = nc.get('v',    idt, slice(None))

    # U/V -> RHO (NaN 보존)
    u=pu.uv2rho_rutgers_safenan(u,"u")
    v=pu.uv2rho_rutgers_safenan(v,"v")
    ubar=pu.uv2rho_rutgers_safenan(ubar,"u")
    vbar=pu.uv2rho_rutgers_safenan(vbar,"v")

    field = ut.ConfigObject(zeta=zeta, ubar=None, vbar=None, temp=temp, salt=salt, u=u, v=v)

# =========================
# 고속/일관 보간 캐시 함수들
# =========================
def build_lni_cache(XD2D, YD2D, Dmask2D, XR2D, YR2D, knn_k=8):
    """
    Donor 바다점만으로 Delaunay 한 번, Target 모든 점의:
    - 삼각형/배리센트릭 가중치(inside)
    - k-NN 인덱스 (outside & 결손 폴백)
    캐시 생성.
    """
    XD2D = np.asarray(XD2D); YD2D = np.asarray(YD2D)
    XR2D = np.asarray(XR2D); YR2D = np.asarray(YR2D)
    src_ok = (np.asarray(Dmask2D) >= 0.5)

    P = np.column_stack((XD2D[src_ok].ravel().astype(np.float64),
                         YD2D[src_ok].ravel().astype(np.float64)))
    tri = Delaunay(P)

    T = np.column_stack((XR2D.ravel().astype(np.float64),
                         YR2D.ravel().astype(np.float64)))
    simp = tri.find_simplex(T)
    inside = (simp >= 0)

    # 배리센트릭 좌표 준비
    trans = tri.transform[simp[inside]]             # (n_in, 3, 2)
    b = np.einsum('nij,nj->ni', trans[:, :2], T[inside] - trans[:, 2])
    w0, w1 = b[:, 0], b[:, 1]
    w2 = 1.0 - w0 - w1
    verts = tri.simplices[simp[inside]]             # (n_in, 3)

    # k-NN 폴백 (outside 및 내부 결손 처리에 공통 사용)
    tree = cKDTree(P)
    _, knn_idx = tree.query(T, k=min(knn_k, len(P)))

    return {
        'src_ok': src_ok,
        'shape': XR2D.shape,
        'inside': inside,
        'verts': verts,
        'w0': w0, 'w1': w1, 'w2': w2,
        'knn_idx': knn_idx,
        'P': P  # 참조용
    }

def apply_lni_cache(values2d, cache):
    """
    캐시를 사용하여 2D 보간 수행.
    - 내부: 배리센트릭 가중평균 (꼭짓점 NaN 있으면 유효한 가중만 재정규화)
    - 외부/결손: k-NN에서 첫 유효값 선택
    """
    vals = np.asarray(values2d, dtype=np.float64)
    vals = vals[cache['src_ok']].ravel()  # donor 바다점 순서로 1D

    N = cache['inside'].size
    out = np.full(N, np.nan, dtype=np.float64)

    # 내부점 처리
    ins = cache['inside']
    if np.any(ins):
        v0 = vals[cache['verts'][:, 0]]
        v1 = vals[cache['verts'][:, 1]]
        v2 = vals[cache['verts'][:, 2]]

        f0 = np.isfinite(v0)
        f1 = np.isfinite(v1)
        f2 = np.isfinite(v2)

        w0 = cache['w0']; w1 = cache['w1']; w2 = cache['w2']

        # 유효값만 가중, 가중합으로 정규화
        num = (np.where(f0, w0*v0, 0.0) +
               np.where(f1, w1*v1, 0.0) +
               np.where(f2, w2*v2, 0.0))
        den = (np.where(f0, w0, 0.0) +
               np.where(f1, w1, 0.0) +
               np.where(f2, w2, 0.0))

        out[ins] = np.divide(num, den, out=np.full_like(num, np.nan), where=den>0)

    # 폴백 (외부 + 내부에서 den==0 등 결손)
    need_fb = ~np.isfinite(out)
    if np.any(need_fb):
        cand_idx = cache['knn_idx'][need_fb]          # (M, k)
        cand_vals = vals[cand_idx]                    # (M, k)
        cand_fin  = np.isfinite(cand_vals)
        # 첫 유효 후보의 위치
        if cand_fin.ndim == 1:
            # k=1인 특수 케이스
            chosen = np.where(cand_fin, cand_vals, np.nan)
        else:
            first = np.argmax(cand_fin, axis=1)
            has  = cand_fin[np.arange(cand_fin.shape[0]), first]
            chosen = np.full(cand_fin.shape[0], np.nan, dtype=np.float64)
            chosen[has] = cand_vals[np.arange(cand_vals.shape[0]), first][has]
        out[need_fb] = chosen

    return out.reshape(cache['shape'])

# --- 캐시 생성 (한 번만) ---
cache_rho = build_lni_cache(lon_src, lat_src, Dmask, lon_dst, lat_dst)

# --- Horizontal interp (2D: Rmask 적용, 3D: Rmask 적용하지 않음) ---
for name in ["zeta","temp","salt","u","v"]:
    print("Horizontal intrp: "+name)
    if not hasattr(field, name):
        continue
    arr = getattr(field, name)
    if arr is None:
        continue

    if arr.ndim == 2:
        out = apply_lni_cache(arr.astype(np.float64, copy=False), cache_rho)
        # 2D는 여기서 Rmask 적용 가능
        out = np.where(Rmask==1, out, 0.0)
        setattr(field, name, out.astype(np.float64, copy=False))

    elif arr.ndim == 3:
        nz = arr.shape[0]
        out = np.empty((nz,) + lon_dst.shape, dtype=np.float64)
        for k in range(nz):
            out[k] = apply_lni_cache(arr[k].astype(np.float64, copy=False), cache_rho)
        # 3D는 아직 Rmask 적용하지 않음 (수직보간 뒤 한 번에)
        setattr(field, name, out)

# --- zr (donor z_r) -> dst 수평보간 (Rmask 적용하지 않음) ---
nz = zr_src.shape[0]
zr_src_remapped = np.empty((nz,) + lon_dst.shape, dtype=np.float64)
for k in range(nz):
    zr_src_remapped[k] = apply_lni_cache(zr_src[k].astype(np.float64, copy=False), cache_rho)

# --- Vertical interpolation to zr_dst ---
field2 = ut.ConfigObject()
for name in ["temp","salt","u","v"]:
    print('Vertical intrp: '+name)
    if hasattr(field, name):
        var = getattr(field, name)
        if var is None or var.ndim != 3:
            continue
        # MATLAB 경계복제와 유사한 padding 외삽 사용(너가 쓰던 설정 유지)
        var_z = pu.vertical_interp_to_ZR(
            zr_src_remapped.astype(np.float64, copy=False),
            var.astype(np.float64, copy=False),
            zr_dst.astype(np.float64, copy=False),
            n_jobs=-1, dedup="mean", extrap_mode="padding",
            zsur=0.0, zbot=-5500.0
        )
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
angle_par_on_dst = apply_lni_cache(angle_src.astype(np.float64, copy=False), cache_rho)

# --- Rotate vectors at rho, then stagger to U/V, then mask_u/v ---
u, v = field2.u, field2.v
Urho, Vrho = ut.rotate_vector_euler(u, v, angle_par_on_dst, to_geo=True)
Ugeo, Vgeo = ut.rotate_vector_euler(Urho, Vrho, angle_dst, to_geo=False)

mask_u = ncG_dst["mask_u"][:]
mask_v = ncG_dst["mask_v"][:]

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

# NaN-aware ρ→U/V 평균 (필요시 ut.rho2uv 교체)
def rho2uv_nanaware(F, pos='u'):
    is2d = (F.ndim == 2)
    if is2d:
        F = F[np.newaxis, ...]
    if pos == 'u':
        a = F[..., :, :-1]
        b = F[..., :,  1:]
    elif pos == 'v':
        a = F[..., :-1, :]
        b = F[...,  1:, :]
    else:
        raise ValueError("pos must be 'u' or 'v'")
    fa = np.isfinite(a); fb = np.isfinite(b)
    out = np.where(fa & fb, 0.5*(a+b),
          np.where(fa, a, np.where(fb, b, np.nan)))
    return np.squeeze(out)

u_on_u = rho2uv_nanaware(Ugeo, pos="u")
v_on_v = rho2uv_nanaware(Vgeo, pos="v")

# 마지막에만 mask_u/v 적용
u_on_u = np.where(mask_u == 1, u_on_u, 0.0)
v_on_v = np.where(mask_v == 1, v_on_v, 0.0)

field2.u = u_on_u.astype(np.float64, copy=False)
field2.v = v_on_v.astype(np.float64, copy=False)

# --- Barotropic from 3D u,v and Hz ---
Hz_dst = np.diff(zw_dst, axis=0).astype(np.float64, copy=False)
ubar, vbar = uv_barotropic_from_3d(field2.u, field2.v, Hz_dst, mask_u=mask_u, mask_v=mask_v)
field2.ubar = ubar.astype(np.float64, copy=False)
field2.vbar = vbar.astype(np.float64, copy=False)

# --- Create initial NetCDF and write ---
status = cn.create_ini__(cfg, grd, 0, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
if status:
    raise RuntimeError(f"Failed creating file {cfg.ininame}")

with Dataset(cfg.ininame, mode='a') as nc:
    nc['zeta'][0] = field2.zeta.astype(nc['zeta'].dtype, copy=False)
    nc['temp'][0] = field2.temp.astype(nc['temp'].dtype, copy=False)
    nc['salt'][0] = field2.salt.astype(nc['salt'].dtype, copy=False)
    nc['u'][0]    = field2.u.astype(nc['u'].dtype, copy=False)
    nc['v'][0]    = field2.v.astype(nc['v'].dtype, copy=False)
    nc['ubar'][0] = field2.ubar.astype(nc['ubar'].dtype, copy=False)
    nc['vbar'][0] = field2.vbar.astype(nc['vbar'].dtype, copy=False)
