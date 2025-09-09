# --- [00] Imports and path setup ---
import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
#import create_I as cn
import utils as ut
import create_I as cn
import post_utils as pu
from scipy.interpolate import LinearNDInterpolator, griddata

pthG_src="/home/shjo/data/roms_inputs/grd/lte/NWP4_grd_3_10m_LP.nc"
pthG_dst="/home/shjo/data/roms_inputs/grd/lte/NWP12_grd_NWP4.nc"
pthI_src="/home/shjo/data/lab/NWP4_ini_DA_8401.nc"
pthI_dst="/home/shjo/data/lab/NWP4_ini_DA_8401.nc"

weight_file="./test_wght.nc"

cfg  = ut.parse_config("./config_all.yaml")
grd  = ut.load_roms_grid(cfg.grdname)

# Read grid
ncG_src=Dataset(pthG_src)
ncG_dst=Dataset(pthG_dst)
ncI_src=Dataset(pthI_src)
ncI_dst=Dataset(pthI_dst)

lon_src,lat_src=ncG_src['lon_rho'][:], ncG_src['lat_rho'][:]
lon_dst,lat_dst=ncG_dst['lon_rho'][:], ncG_dst['lat_rho'][:]

topo_src = ncG_src["h"][:]
topo_dst = ncG_dst["h"][:]
mask_dst = ncG_dst["mask_rho"][:]

angle_src = ncG_src["angle"][:]
angle_dst = ncG_dst["angle"][:]

vtrs=float(ncI_src["Vtransform"][:].data)
vstr=float(ncI_src["Vstretching"][:].data)
theta_s=float(ncI_src["theta_s"][:].data)
theta_b=float(ncI_src["theta_b"][:].data)
theta_b=0.1
tcline=float(ncI_src["Tcline"][:].data)
nlayer=len(ncI_src["sc_r"][:])

zr_src=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,nlayer,1,topo_src,np.zeros_like(topo_src))
zr_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,36,1,topo_dst,np.zeros_like(topo_dst))

zw_dst=ut.zlevs(vtrs,vstr,theta_s,theta_b,tcline,36,5,topo_dst,np.zeros_like(topo_dst))

# --- crop domain ---
lon_crop, lat_crop, idx, idy = ut.crop_to_model_domain(lat_dst, lon_dst, lat_src, lon_src)

# --- calc weight ---
status = ut.build_bilinear_regridder(lon_src, lat_src, lon_dst, lat_dst, weight_file, reuse=False)
if status:
    raise RuntimeError(f"Failed to generate remap weights: {weight_file}")

# --- Load & calc increments ---
factor=1
zeta_diff=(ncI_src["zeta"][0]-ncI_src["zeta"][1])*factor
temp_diff=(ncI_src["temp"][0]-ncI_src["temp"][1])*factor
salt_diff=(ncI_src["salt"][0]-ncI_src["salt"][1])*factor
u_diff=(ncI_src["u"][0]-ncI_src["u"][1])*factor
v_diff=(ncI_src["v"][0]-ncI_src["v"][1])*factor
ubar_diff=(ncI_src["ubar"][0]-ncI_src["ubar"][1])*factor
vbar_diff=(ncI_src["vbar"][0]-ncI_src["vbar"][1])*factor

u_diff=pu.uv2rho_rutgers_safenan(u_diff,"u")
v_diff=pu.uv2rho_rutgers_safenan(v_diff,"v")
ubar_diff=pu.uv2rho_rutgers_safenan(ubar_diff,"u")
vbar_diff=pu.uv2rho_rutgers_safenan(vbar_diff,"v")

def interp_angle_like_matlab(angle_src, lon_src, lat_src, lon_dst, lat_dst):
    """
    MATLAB roms2roms의 angle 처리와 동일한 효과:
    - angle을 '선형 보간'으로 dst 격자에 내삽
    - convex hull 밖/결측은 'nearest'로 보완
    """
    # 1) 선형 보간 (MATLAB scatteredInterpolant linear와 등가)
    valid = ~np.isnan(angle_src)
    if valid.sum() == 0:
        # 소스 전체 NaN이면 그대로 반환
        return np.full_like(lon_dst, np.nan, dtype=angle_src.dtype)

    interp_lin = LinearNDInterpolator(
        (lon_src[valid], lat_src[valid]),
        angle_src[valid],
        fill_value=np.nan
    )
    ang_dst = interp_lin(lon_dst, lat_dst)

    # 2) 남은 NaN만 최근접으로 보완 (MATLAB에서 외곽/구멍 fallback)
    nan_mask = np.isnan(ang_dst)
    if np.any(nan_mask):
        points  = np.column_stack((lon_src[valid], lat_src[valid]))
        values  = angle_src[valid]
        targets = (lon_dst[nan_mask], lat_dst[nan_mask])
        ang_dst[nan_mask] = griddata(points, values, targets, method='nearest')

    # (옵션) dtype 유지
    return ang_dst.astype(angle_src.dtype, copy=False)

def hz_to_u_v(Hz: np.ndarray):
    """
    Hz @RHO (nz, ny, nx) -> 층두께를 C-grid U/V 위치로 수평 평균
    dzu: (nz, ny, nx-1),  dzv: (nz, ny-1, nx)
    """
    # U-방향 평균 (x-방향 이웃 평균)
    dzu = 0.5 * (Hz[:, :, :-1] + Hz[:, :, 1:])
    # V-방향 평균 (y-방향 이웃 평균)
    dzv = 0.5 * (Hz[:, :-1, :] + Hz[:, 1:, :])
    return dzu, dzv

def uv_barotropic_from_3d(u: np.ndarray,
                          v: np.ndarray,
                          Hz: np.ndarray,
                          mask_u: np.ndarray | None = None,
                          mask_v: np.ndarray | None = None,
                          eps: float = 1e-12):
    """
    MATLAB uv_barotropic (full grid)과 동일한 수식으로 barotropic 계산.

    입력
      u : (nz, ny, nx-1)  @U-point
      v : (nz, ny-1, nx)  @V-point
      Hz: (nz, ny, nx)    @RHO-point (층두께; 타깃 격자)
      mask_u, mask_v: (ny, nx-1)/(ny-1, nx)  (선택) 최종 마스킹용

    출력
      ubar: (ny, nx-1)
      vbar: (ny-1, nx)
    """
    # RHO->U/V 위치로 층두께 변환 (MATLAB의 Duk/Dvk 정의와 동일)
    dzu, dzv = hz_to_u_v(Hz)

    # 분모(수심 적분 두께)
    Du = np.sum(dzu, axis=0)            # (ny, nx-1)
    Dv = np.sum(dzv, axis=0)            # (ny-1, nx)

    # 분자(가중합)
    num_u = np.sum(u * dzu, axis=0)     # (ny, nx-1)
    num_v = np.sum(v * dzv, axis=0)     # (ny-1, nx)

    # 가중평균 (MATLAB: ubar=ubar./Du; vbar=... 와 동일)
    ubar = np.where(Du > eps, num_u / Du, 0.0)
    vbar = np.where(Dv > eps, num_v / Dv, 0.0)

    # (선택) 마스크 적용: MATLAB은 마지막에 경계 바깥이 0이 되도록 처리하는 관행
    if mask_u is not None:
        ubar = np.where(mask_u == 1, ubar, 0.0)
    if mask_v is not None:
        vbar = np.where(mask_v == 1, vbar, 0.0)

    return ubar, vbar

# --- define field ---
field = ut.ConfigObject(zeta=zeta_diff, temp=temp_diff, salt=salt_diff,\
        u=u_diff,v=v_diff,ubar=ubar_diff,vbar=vbar_diff)

# --- Load and apply remap weights to all fields ---
with Dataset(weight_file) as nc:
    row = nc.variables["row"][:] - 1
    col = nc.variables["col"][:] - 1
    S   = nc.variables["S"][:]
for varname in vars(field):
    var_src = getattr(field, varname)
    remapped = ut.remap_variable(var_src, row, col, S, lon_dst.shape, method="coo")
    setattr(field, varname, remapped)

# --- Horizontal flood (all fields) ---
for var in vars(field):
    val = getattr(field, var)
    val_flooded = ut.flood_horizontal(val, lon_dst, lat_dst, method="griddata")
    setattr(field, var, val_flooded)
#    done("flood_h", time.time()-t0)

# --- interp zr ---
zr_src_remapped = ut.remap_variable(zr_src, row, col, S, lon_dst.shape, method="coo")
zr_src_flooded = ut.flood_horizontal(zr_src_remapped, lon_dst, lat_dst, method="griddata")
"""
# --- Mask land to 0 ---
for varname in vars(field):
    var = getattr(field, varname)
    if var.ndim==2:
        var[mask_dst == 0] = 0.0
    else:
        var[:, mask_dst == 0] = 0.0
    setattr(field, varname, var)
"""
# --- vertical interpolation ---
for varname in vars(field):
    var = getattr(field,varname)
    if var.ndim==2:
        continue
    print(varname)
    var = getattr(field, varname)
    var_zinterp = pu.vertical_interp_to_ZR(zr_src_flooded, var, zr_dst,
                                         n_jobs=-1,
                                         dedup="mean",
                                         extrap_mode="padding")
    setattr(field, varname, var_zinterp)

# angle_src -> T.parent_angle (dst로 선형보간)
angle_par_on_dst = interp_angle_like_matlab(angle_src, lon_src, lat_src, lon_dst, lat_dst)

# --- rotate vector ---
# 1. parent grid -> geo (rho-point 회전)
u,v=getattr(field,"u"), getattr(field,"v")
Urho, Vrho = ut.rotate_vector_euler(u, v, angle_par_on_dst, to_geo=True)

# 2. geo -> target grid (rho-point 회전)
Ugeo, Vgeo = ut.rotate_vector_euler(Urho, Vrho, angle_dst, to_geo=False)

mask_u = ncG_dst["mask_u"][:]
mask_v = ncG_dst["mask_v"][:]

u_on_u = ut.rho2uv(Ugeo, pos="u")
v_on_v = ut.rho2uv(Vgeo, pos="v")
u_on_u = np.where(mask_u == 1, u_on_u, 0.0)
v_on_v = np.where(mask_v == 1, v_on_v, 0.0)
setattr(field, "u", u_on_u)
setattr(field, "v", v_on_v)

print(zw_dst.shape)

Hz_dst = np.diff(zw_dst,axis=0)

print(Hz_dst.shape)
ubar, vbar = uv_barotropic_from_3d(field.u, field.v, Hz_dst, mask_u=mask_u, mask_v=mask_v)

setattr(field, "ubar", ubar)
setattr(field, "vbar", vbar)


# [03] Create initial NetCDF file
status = cn.create_ini__(cfg, grd, 0, ncFormat=cfg.ncformat, bio_model=cfg.bio_model_type)
if status:
    raise RuntimeError(f"Failed creating file {cfg.ininame}")


# [13] Write all variables to ini.nc
with Dataset(cfg.ininame, mode='a') as nc:
    for var in ['zeta', 'temp', 'salt', 'u', 'v', 'ubar', 'vbar']:
        nc[var][0] = getattr(field, var)









