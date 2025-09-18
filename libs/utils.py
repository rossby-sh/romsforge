import numpy as np
import datetime as dt
from netCDF4 import Dataset, date2num
import yaml

class ConfigObject:
    def __init__(self, **entries):
        for key, value in entries.items():
            if isinstance(value, dict):
                value = ConfigObject(**value)
            self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self, indent=0):
        lines = []
        pad = '  ' * indent
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigObject):
                lines.append(f"{pad}{key}:")
                lines.append(value.__repr__(indent + 1))
            else:
                lines.append(f"{pad}{key}: {value}")
        return '\n'.join(lines)

    def to_dict(self):
        out = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigObject):
                out[key] = value.to_dict()
            else:
                out[key] = value
        return out

class MaskedNetCDF:
    def __init__(self, nc, atol=1e-6):
        self.nc = nc
        self.atol = atol

    def get(self, name, *slices):
        var = self.nc[name][*slices]
        var = var.astype(float, copy=False)

        # MaskedArray 처리
        if isinstance(var, np.ma.MaskedArray) and np.any(var.mask):
            var = var.filled(np.nan)

        # _FillValue 처리
        if hasattr(self.nc[name], '_FillValue'):
            fv = self.nc[name]._FillValue
            mask = np.isclose(var, fv, atol=self.atol)
            var[mask] = np.nan

        return var


def parse_config(path="config.yaml") -> ConfigObject:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return ConfigObject(**raw)

def compute_relative_time(time, input_unit, ref_unit):
    """
    시간 데이터를 ref_unit 기준 상대 시간으로 변환

    Parameters:
        time        : 원본 시간 배열 (np.ndarray)
        input_unit  : OGCM 시간 단위 (ex: 'hours since ...', 'seconds since ...')
        ref_unit    : ROMS 기준 시간 단위 (ex: 'days since ...')

    Returns:
        np.ndarray : ref_unit 기준 상대 시간
    """
    anchor = dt.datetime(2025, 1, 1)

    # 기준 시간 계산
    t_in  = date2num(anchor, input_unit)
    t_ref = date2num(anchor, ref_unit)

    # 둘 다 days 기준으로 변환
    def to_days(value, unit):
        if "seconds" in unit:
            return value / 86400.0
        elif "hours" in unit:
            return value / 24.0
        elif "days" in unit:
            return value
        else:
            raise ValueError(f"Unknown time unit: {unit}")

    t_in_days = to_days(t_in, input_unit)
    t_ref_days = to_days(t_ref, ref_unit)
    offset = t_in_days - t_ref_days

    # time도 days로 변환
    time_in_days = to_days(time, input_unit)
    total_days = time_in_days - offset

    # 출력 단위 맞추기
    if "seconds" in ref_unit:
        return total_days * 86400.0
    elif "hours" in ref_unit:
        return total_days * 24.0
    elif "days" in ref_unit:
        return total_days
    else:
        raise ValueError(f"Unknown ref_unit: {ref_unit}")


def crop_to_model_domain(lat_src, lon_src, lat_dst, lon_dst):
    """
    lat/lon이 1D (rectilinear) 또는 2D (curvilinear, WRF)인 경우 모두 지원.

    Returns:
        idx, idy: x, y 방향 인덱스 범위
        lon_crop, lat_crop: 잘린 source grid (1D 또는 2D 형태 그대로 유지)
    """

    if lat_src.ndim == 1 and lon_src.ndim == 1:
        lon_coarse = np.where((lon_src >= np.min(lon_dst)) & (lon_src <= np.max(lon_dst)))[0]
        lat_coarse = np.where((lat_src >= np.min(lat_dst)) & (lat_src <= np.max(lat_dst)))[0]

        lon_ref = lon_src[lon_coarse]
        lat_ref = lat_src[lat_coarse]

        lon_res = max(
            np.max(np.abs(np.diff(lon_ref))),
            np.max(np.abs(np.diff(lon_dst[0, :])))
        )
        lat_res = max(
            np.max(np.abs(np.diff(lat_ref))),
            np.max(np.abs(np.diff(lat_dst[:, 0])))
        )

        lon_min = min(np.min(lon_ref), np.min(lon_dst) - lon_res)
        lon_max = max(np.max(lon_ref), np.max(lon_dst) + lon_res)
        lat_min = min(np.min(lat_ref), np.min(lat_dst) - lat_res)
        lat_max = max(np.max(lat_ref), np.max(lat_dst) + lat_res)

        idx = np.where((lon_src >= lon_min) & (lon_src <= lon_max))[0]
        idy = np.where((lat_src >= lat_min) & (lat_src <= lat_max))[0]

        lon_crop = lon_src[idx]
        lat_crop = lat_src[idy]

        lon_crop, lat_crop = np.meshgrid(lon_crop, lat_crop)

        return  lon_crop, lat_crop, idx, idy

    elif lat_src.ndim == 2 and lon_src.ndim == 2:
        lat_min, lat_max = lat_dst.min(), lat_dst.max()
        lon_min, lon_max = lon_dst.min(), lon_dst.max()

        lat_res = np.max(np.abs(np.diff(lat_src[:, 0])))
        lon_res = np.max(np.abs(np.diff(lon_src[0, :])))

        lat_min -= lat_res
        lat_max += lat_res
        lon_min -= lon_res
        lon_max += lon_res

        mask = (
            (lat_src >= lat_min) & (lat_src <= lat_max) &
            (lon_src >= lon_min) & (lon_src <= lon_max)
        )

        j_idx, i_idx = np.where(mask)
        j_min, j_max = j_idx.min(), j_idx.max() + 1
        i_min, i_max = i_idx.min(), i_idx.max() + 1

        idy = np.arange(j_min, j_max)
        idx = np.arange(i_min, i_max)

        lon_crop = lon_src[np.ix_(idy, idx)]
        lat_crop = lat_src[np.ix_(idy, idx)]

        return lon_crop, lat_crop, idx, idy


def load_roms_grid(grdname):
    """
    ROMS 그리드 파일에서 필요한 변수들만 추출해 dict 형태로 반환

    Returns:
        dict: {"lon", "lat", "angle", "h", "mask"}
    """
    with Dataset(grdname) as nc:
        lon = nc["lon_rho"][:]
        lat = nc["lat_rho"][:]
        angle = nc["angle"][:]
        h = nc["h"][:]
        mask_rho = nc["mask_rho"][:]
        mask_u   = nc["mask_u"][:]
        mask_v   = nc["mask_v"][:]

    return ConfigObject(
        lon     = lon,
        lat     = lat,
        angle   = angle,
        topo    = h,
        mask_rho= mask_rho,
        mask_u  = mask_u,
        mask_v  = mask_v

    )

def load_ogcm_metadata(ogcm_name,ogcm_var_name):
    """
    OGCM 파일에서 위경도, 깊이, 시간 정보를 읽어와 dict로 반환

    Parameters:
        ogcm_path (str): OGCM NetCDF 파일 경로
        varmap (dict): config.yaml에서 정의한 변수명 매핑

    Returns:
        dict: {"lon", "lat", "depth", "time"}
    """
    with Dataset(ogcm_name) as nc:
        lon = nc[ogcm_var_name["longitude"]][:]
        lat = nc[ogcm_var_name["latitude"]][:]
        depth = nc[ogcm_var_name["depth"]][:]
        time = nc[ogcm_var_name["time"]][:]
        time_unit = nc[ogcm_var_name["time"]].units
    return ConfigObject(
        lon       = lon,
        lat       = lat,
        depth     = depth,
        time      = time,
        time_unit = time_unit
    )

def build_bilinear_regridder(lon_src, lat_src, lon_dst, lat_dst, wghtname, reuse=False):
    """XESMF 래퍼 함수: 감정 없는 weight 생성기"""
    import xarray as xr
    import xesmf as xe
    try:
        ds_src = xr.Dataset({
            "lon": (["y", "x"], lon_src),
            "lat": (["y", "x"], lat_src)
        })
        ds_dst = xr.Dataset({
            "lon": (["y", "x"], lon_dst),
            "lat": (["y", "x"], lat_dst)
        })

        xe.Regridder(ds_src, ds_dst, method="bilinear", filename=wghtname, reuse_weights=reuse)
        print(f"--- [+] Weight file created: {wghtname} ---")
        return 0

    except Exception as e:
        print(f"--- [!ERROR] Failed to create weight file: {e} ---")
        return 1


def depth_average(var3d, depth, mask_nan=True):
    """
    depth 방향으로 평균 내는 함수

    Parameters
    ----------
    var3d : ndarray
        (z, y, x) 형식의 3D 변수 (ex: u, v)
    depth : 1D array
        depthO 또는 z-level (길이 z)
    mask_nan : bool
        NaN 값 제외할지 여부

    Returns
    -------
    avg2d : ndarray
        depth-averaged 2D 필드 (y, x)
    """
    dz = np.gradient(-depth)
    du = np.zeros_like(var3d[0])
    zu = np.zeros_like(var3d[0])

    for n in range(len(depth)):
        layer = var3d[n]
        valid = ~np.isnan(layer) if mask_nan else np.ones_like(layer)
        du += dz[n] * np.where(np.isnan(layer), 0, layer)
        zu += dz[n] * valid
    avg = du / zu
    avg[zu == 0] = np.nan
    return avg


def remap_variable(var_src, row, col, S, dst_shape, method="coo"):
    """
    SCRIP 가중치를 사용한 2D 또는 3D 변수 리매핑 (기본: COO 방식)

    Parameters
    ----------
    var_src : ndarray
        2D 또는 3D 원본 변수 (e.g., zeta, temp 등)
    row, col, S : ndarray
        SCRIP weight 정보
    dst_shape : tuple
        출력할 데이터 lon_rho shape ( (y, x))
    method : str
        리매핑 방식 (기본: "coo", 추후 "crs" 등 추가 가능)

    Returns
    -------
    var_dst : ndarray
        리매핑된 출력 변수
    """
    if method != "coo":
        raise NotImplementedError(f"remap method '{method}' is not implemented yet.")

    if var_src.ndim == 2:
        ny_src, nx_src = var_src.shape
        ydim, xdim = dst_shape
        src_flat = var_src.flatten()
        dst_flat = np.zeros(row.max() + 1)

        np.add.at(dst_flat, row, S * src_flat[col])
        var_dst = dst_flat.reshape((ydim, xdim))

    elif var_src.ndim == 3:
        nz, ny_src, nx_src = var_src.shape
        ydim, xdim = dst_shape
        src_flat = var_src.reshape(nz, ny_src * nx_src)
        var_dst = np.zeros((nz, ydim, xdim))

        for k in range(nz):
            dst_flat = np.zeros(row.max() + 1)
            np.add.at(dst_flat, row, S * src_flat[k, col])
            var_dst[k] = dst_flat.reshape((ydim, xdim))

    else:
        raise ValueError("Only 2D or 3D input arrays are supported.")

    return var_dst

def rotate_vector_euler(u, v, angle, to_geo=True):
    is2d = (u.ndim == 2)

    if is2d:
        u = u[np.newaxis, ...]
        v = v[np.newaxis, ...]

    nz, ny, nx = u.shape
    angle3d = np.broadcast_to(angle, (nz, ny, nx))

    uv_complex = u + 1j * v

    if to_geo:
        uv_rotated = uv_complex * np.exp(-1j * angle3d)
    else:
        uv_rotated = uv_complex * np.exp(+1j * angle3d)

    return np.squeeze(uv_rotated.real), np.squeeze(uv_rotated.imag)

def rho2uv(field, pos='u'):
    """
    Convert rho-point data to u-point or v-point.

    Parameters:
        field: np.ndarray (2D or 3D) — shape (ny, nx) or (nz, ny, nx)
        pos: 'u' or 'v' — direction to convert

    Returns:
        np.ndarray of shape [..., ny, nx-1] or [..., ny-1, nx]
    """
    # Ensure at least 3D shape: (nz, ny, nx)
    if field.ndim == 2:
        field = field[np.newaxis, ...]

    if pos == 'u':
        result = 0.5 * (field[..., :, :-1] + field[..., :, 1:])
    elif pos == 'v':
        result = 0.5 * (field[..., :-1, :] + field[..., 1:, :])
    else:
        raise ValueError("pos must be 'u' or 'v'")

    return np.squeeze(result)


def conserve_and_recompute_barotropic(u, v, ubar, vbar, dzu, dzv):
    """Apply volume conservation correction and recompute barotropic velocities."""
    u_corr = u - np.sum(u * dzu, axis=0) / np.sum(dzu, axis=0) + ubar
    v_corr = v - np.sum(v * dzv, axis=0) / np.sum(dzv, axis=0) + vbar

    ubar_new = np.sum(u_corr * dzu, axis=0) / np.sum(dzu, axis=0)
    vbar_new = np.sum(v_corr * dzv, axis=0) / np.sum(dzv, axis=0)

    return u_corr, v_corr, ubar_new, vbar_new


def stretching(Vstretching,theta_s,theta_b,Layer_N,kgrid=1):
    
    # Vstretching=MyVar['Vstretching']
    # theta_s=MyVar['Theta_s']
    # theta_b=MyVar['Theta_b']
    # N=MyVar['Layer_N']    
    Np=Layer_N+1
    if Vstretching==1:
        ds=1/Layer_N
        if kgrid==1:
            nlev=Np
            lev=np.arange(Layer_N+1)
            s=(lev-Layer_N)*ds
        else:
            Nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-.5
            s=(lev-Layer_N)*ds
        
        if theta_s>0:
            Ptheta=np.sinh(theta_s*s)/np.sinh(theta_s)
            Rtheta=np.tanh(theta_s*(s+0.5))/(2.0*np.tanh(0.5*theta_s))-0.5
            C=(1.0-theta_b)*Ptheta+theta_b*Rtheta
        else:
            C=s 
        
    elif Vstretching==2:
        
        alfa=1.0
        beta=1.0
        ds=1.0/Layer_N
        if kgrid==1:
            Nlev=Np
            lev=np.arange(Layer_N+1)
            s=(lev-Layer_N)*ds
        else:
            Nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-.5
            s=(lev-Layer_N)*ds
        
        if theta_s>0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
            if theta_b>0:
                Cbot=-1.0+np.sinh(theta_b*(s+1.0))/np.sinh(theta_b)
                weigth=(s+1.0)**alfa*(1.0+(alfa/beta)*(1.0-(s+1)**beta))
                C=weigth*Csur+(1.0-weigth)*Cbot
            else:
                C=Csur
    elif Vstretching==4:
        ds=1.0/Layer_N
        if kgrid==1:
            Nlev=Np
            lev=np.arange(Layer_N+1)
            s=(lev-Layer_N)*ds
        else:
            nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-0.5
            s=(lev-Layer_N)*ds
        if theta_s>0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            Csur=-s**2
        
        if theta_b>0:
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
        
        
    elif Vstretching==5:
        if kgrid==1:
            nlev=Np
            lev=np.arange(Layer_N+1)
            s=-(lev*lev-2.0*lev*Layer_N+lev+Layer_N*Layer_N-Layer_N)/(Layer_N*Layer_N-Layer_N)-\
                0.01*(lev*lev-lev*Layer_N)/(1.0-Layer_N)
            s[0]=-1.0
        else:
            Nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-0.5
            s=-(lev*lev-2.0*lev*Layer_N+lev+Layer_N*Layer_N-Layer_N)/(Layer_N*Layer_N-Layer_N)-\
                0.01*(lev*lev-lev*Layer_N)/(1.0-Layer_N)
        if theta_s>0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            Csur=-s**2
        if theta_b>0:
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
    return s,C

from numba import jit 
@jit(nopython=True,cache=True) 
def ztosigma_numba(var, z, depth):
    Ns, Mp, Lp = z.shape
    Nz = len(depth)
    vnew = np.zeros((Ns, Mp, Lp))

    for ks in range(Ns):
        sigmalev = z[ks]
        thezlevs = np.zeros((Mp, Lp), dtype=np.int32)

        for j in range(Mp):
            for i in range(Lp):
                for kz in range(Nz):
                    if sigmalev[j, i] > depth[kz]:
                        thezlevs[j, i] += 1

        min_lev = thezlevs.min()
        max_lev = thezlevs.max()
        if max_lev >= Nz or min_lev <= 0:
            break

        for j in range(Mp):
            for i in range(Lp):
                kz = thezlevs[j, i]
                z1 = depth[kz - 1]
                z2 = depth[kz]
                v1 = var[kz - 1, j, i]
                v2 = var[kz, j, i]
                sigma_z = sigmalev[j, i]
                vnew[ks, j, i] = ((v1 - v2) * sigma_z + v2 * z1 - v1 * z2) / (z1 - z2)

    return vnew

from numba import njit, prange

@njit(parallel=True)
def ztosigma_numba_parallel(var, z, depth):
    Ns, Mp, Lp = z.shape
    Nz = len(depth)
    vnew = np.zeros((Ns, Mp, Lp))

    for ks in range(Ns):
        sigmalev = z[ks]
        thezlevs = np.zeros((Mp, Lp), dtype=np.int32)

        for j in range(Mp):
            for i in range(Lp):
                for kz in range(Nz):
                    if sigmalev[j, i] > depth[kz]:
                        thezlevs[j, i] += 1

        min_lev = thezlevs.min()
        max_lev = thezlevs.max()
        if max_lev >= Nz or min_lev <= 0:
            break

        for j in range(Mp):
            for i in range(Lp):
                kz = thezlevs[j, i]
                z1 = depth[kz - 1]
                z2 = depth[kz]
                v1 = var[kz - 1, j, i]
                v2 = var[kz, j, i]
                sigma_z = sigmalev[j, i]
                vnew[ks, j, i] = ((v1 - v2) * sigma_z + v2 * z1 - v1 * z2) / (z1 - z2)

    return vnew

@jit(nopython=True)
def ztosigma_1d_numba(var, z, depth):
    Ns, Lp = z.shape
    Nz = len(depth)
    vnew = np.zeros((Ns, Lp))
    for ks in range(Ns):
        sigmalev = z[ks, :]
        thezlevs = np.zeros(Lp, dtype=np.int32)
        for kz in range(Nz):
            for j in range(Lp):
                if sigmalev[j] > depth[kz]:
                    thezlevs[j] += 1
        for j in range(Lp):
            if thezlevs[j] >= Nz or thezlevs[j] <= 0:
                continue
            k1 = thezlevs[j] - 1
            k2 = thezlevs[j]
            z1 = depth[k1]
            z2 = depth[k2]
            tmp_var = var[:, j]
            v1 = tmp_var[k1]
            v2 = tmp_var[k2]
            vnew[ks, j] = ((v1 - v2) * sigmalev[j] + v2 * z1 - v1 * z2) / (z1 - z2)
    return vnew


def ztosigma(var,z,depth):
    Ns,Mp,Lp=z.shape
    Nz=len(depth)
    vnew=np.zeros([Ns,Mp,Lp])
    for ks in range(Ns):
        sigmalev=np.squeeze(z[ks,:,:])
        thezlevs=0*sigmalev
        for kz in range(Nz):
            thezlevs[sigmalev>depth[kz]]=thezlevs[sigmalev>depth[kz]]+1
        if np.max(thezlevs)>=Nz or np.min(thezlevs)<=0:
            print("min sigma level = "+str(np.min(z))+' - min z level = '+\
                  str(np.min(depth)))
            print("max sigma level = "+str(np.max(z))+' - max z level = '+\
                  str(np.max(depth)))            
        thezlevs=thezlevs.astype('int32')
        imat,jmat=np.meshgrid(np.arange(1,Lp+1),np.arange(1,Mp+1))
        pos=Nz*Mp*(imat-1)+Nz*(jmat-1)+thezlevs
        z1,z2=depth[thezlevs-1],depth[thezlevs]
        tmp_var=var.transpose().flatten()
        v1=tmp_var[pos-1].reshape(Mp,Lp)
        v2=tmp_var[pos].reshape(Mp,Lp)
        vnew[ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vnew

def ztosigma_1d(var,z,depth):
    Ns,Lp=z.shape
    Nz=len(depth)
    Nz=var.shape[0]
    vnew=np.zeros([Ns,Lp])
    for ks in range(Ns):
        sigmalev=np.squeeze(z[ks,:])
        thezlevs=0*sigmalev
        for kz in range(Nz):
            thezlevs[sigmalev>depth[kz]]=thezlevs[sigmalev>depth[kz]]+1
        if np.max(thezlevs)>=Nz or np.min(thezlevs)<=0:
            print("min sigma level = "+str(np.min(z))+' - min z level = '+\
                  str(np.min(depth)))
            print("max sigma level = "+str(np.max(z))+' - max z level = '+\
                  str(np.max(depth))) 
                
        thezlevs=thezlevs.astype('int32')
        
        jmat= np.arange(1,Lp+1)
        pos=Nz*(jmat-1)+thezlevs
        

        z1,z2=depth[thezlevs-1],depth[thezlevs]
        tmp_var=var.transpose().flatten()

        v1=tmp_var[pos-1].reshape(Lp)
        v2=tmp_var[pos].reshape(Lp)
        vnew[ks,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vnew

        
def zlevs(Vtransform, Vstretching,theta_s, theta_b, hc, N,igrid, h, zeta):
    # Vstretching=2;Vtransform=2;theta_s=7;theta_b=0.1;hc=200;N=36;igrid=1;h=topo_dst; zeta=np.zeros_like(h)
    from copy import deepcopy
    
    # for get section
    if len(h.shape)!=2:
        h=h.reshape([1,len(h)])
        zeta=zeta.reshape([1,len(zeta)])
    Np=N+1;
    Lp,Mp=h.shape
    L=Lp-1;
    M=Mp-1;
    
    hmin=np.min(h);
    hmax=np.max(h);

    # Compute vertical stretching function, C(k):

    if igrid==5:
        kgrid=1
        s,C=stretching(Vstretching, theta_s, theta_b, N, kgrid)
    else:
        kgrid=0
        s,C=stretching(Vstretching, theta_s, theta_b, N, kgrid)

    #  Average bathymetry and free-surface at requested C-grid type.

    if igrid==1:
        hr=deepcopy(h)
        zetar=deepcopy(zeta)
    elif igrid==2:
        hp=0.25*(h[:L,:M]+h[1:Lp,:M]+h[:L,1:Mp]+h[1:Lp,1:Mp])
        zetap=0.25*(zeta[:L,:M]+zeta[1:Lp,:M]+zeta[:L,1:Mp]+zeta[1:Lp,1:Mp])
    elif igrid==3:
        hu=0.5*(h[:L,:Mp]+h[1:Lp,:Mp])
        zetau=0.5*(zeta[:L,:Mp]+zeta[1:Lp,:Mp])
    elif igrid==4:
        hv=0.5*(h[:Lp,:M]+h[:Lp,1:Mp])
        zetav=0.5*(zeta[:Lp,:M]+zeta[:Lp,1:Mp])        
    elif igrid==5:
        hr=deepcopy(h)
        zetar=deepcopy(zeta)


    # Compute depths (m) at requested C-grid location.
    if Vtransform==1:
        if igrid==1:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hr
                z[:,:,k]=z0+zetar*(1+z0/hr)
        elif igrid==2:
            z=np.zeros([zetap.shape[0],zetap.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hp
                z[:,:,k]=z0+zetap*(1+z0/hp)
        elif igrid==3:
            z=np.zeros([zetau.shape[0],zetau.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hu
                z[:,:,k]=z0+zetau*(1+z0/hu)
        elif igrid==4:
            z=np.zeros([zetav.shape[0],zetav.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hv
                z[:,:,k]=z0+zetav*(1+z0/hv)
        elif igrid==5:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],Np])
            z[:,:,0]=-hr
            for k in range(1,Np):
                z0=(s[k]-C[k])*hc+C[k]*hr
                z[:,:,k]=z0+zetar*(1+z0/hr)
    
    elif Vtransform==2:
        if igrid==1:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hr)/(hc+hr)
                z[:,:,k]=zetar+(zeta+hr)*z0
        elif igrid==2:
            z=np.zeros([zetap.shape[0],zetap.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hp)/(hc+hp)
                z[:,:,k]=zetap+(zetap+hp)*z0
        elif igrid==3:
            z=np.zeros([zetau.shape[0],zetau.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hu)/(hc+hu)
                z[:,:,k]=zetau+(zetau+hu)*z0
        elif igrid==4:
            z=np.zeros([zetav.shape[0],zetav.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hv)/(hc+hv)
                z[:,:,k]=zetav+(zetav+hv)*z0
        elif igrid==5:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],Np])
            for k in range(0,Np):
               z0=(hc*s[k]+C[k]*hr)/(hc+hr)
               z[:,:,k]=zetar+(zetar+hr)*z0
    z=np.squeeze(np.transpose(z,[2,0,1]))
    
    return z


from scipy.ndimage import distance_transform_edt
from scipy.interpolate import LinearNDInterpolator, griddata

def flood_horizontal(var: np.ndarray, lon2d: np.ndarray, lat2d: np.ndarray, method: str = 'edt') -> np.ndarray:
    """
    Flood missing (NaN) values in 2D or 3D field using specified horizontal method.
    
    Parameters:
    - var: 2D or 3D ndarray (e.g., temp, salt, u, v)
    - lon2d, lat2d: 2D horizontal grid
    - method: 'edt' or 'linearND'
    
    Returns:
    - var_filled: array with NaNs filled
    """
    if var.ndim == 2:
        return _flood_2d(var, lon2d, lat2d, method)

    elif var.ndim == 3:
        Nz = var.shape[0]
        filled = np.empty_like(var)
        for k in range(Nz):
            filled[k] = _flood_2d(var[k], lon2d, lat2d, method)
        return filled

    else:
        raise ValueError(f"Unsupported dimension: {var.ndim}")


def _flood_2d(var2d: np.ndarray, lon2d: np.ndarray, lat2d: np.ndarray, method: str) -> np.ndarray:
    nan_mask = np.isnan(var2d)
    if not np.any(nan_mask):
        return var2d

    var_filled = var2d.copy()
    valid = ~nan_mask

    if method == 'edt':
        _, indices = distance_transform_edt(nan_mask, return_indices=True)
        var_filled[nan_mask] = var2d[tuple(indices)][nan_mask]

    elif method == 'linearND':
        interp = LinearNDInterpolator(
            (lon2d[valid], lat2d[valid]),
            var2d[valid],
            fill_value=0.0  # fallback
        )
        var_filled[nan_mask] = interp(lon2d[nan_mask], lat2d[nan_mask])

    elif method == 'griddata':
        points = np.column_stack((lon2d[valid], lat2d[valid]))
        values = var2d[valid].flatten()
        all_points = np.column_stack((lon2d.ravel(), lat2d.ravel()))
        interpolated = griddata(points, values, all_points, method='nearest', fill_value=0.0)
        interpolated_2d = interpolated.reshape(var2d.shape)
        var_filled[nan_mask] = interpolated_2d[nan_mask]

    else:
        raise ValueError(f"Unknown method: {method}")

    return var_filled

from scipy.spatial import cKDTree

def _flood_2d_new(var2d, lon2d, lat2d, method):
    var2d = np.asarray(var2d)
    nan_mask = np.isnan(var2d)
    if not np.any(nan_mask):
        return var2d

    var_filled = var2d.copy()
    valid = ~nan_mask
    if valid.sum() == 0:
        return var2d  # all-NaN guard

    if method == 'edt':
        _, indices = distance_transform_edt(nan_mask, return_indices=True)
        var_filled[nan_mask] = var2d[tuple(indices)][nan_mask]

    elif method == 'linearND':
        # 1단계: 선형 보간(볼록껍질 밖은 NaN)
        interp = LinearNDInterpolator(
            (lon2d[valid], lat2d[valid]),
            var2d[valid],
            fill_value=np.nan
        )
        vals = interp(lon2d, lat2d)
        # 2단계: 남은 NaN만 최근접으로 채움
        nan2 = np.isnan(vals)
        if np.any(nan2):
            pts_valid = np.column_stack((lon2d[valid], lat2d[valid]))
            pts_nan   = np.column_stack((lon2d[nan2],  lat2d[nan2]))
            tree = cKDTree(pts_valid)
            _, idx = tree.query(pts_nan, k=1)
            vals[nan2] = var2d[valid][idx]
        var_filled[nan_mask] = vals[nan_mask]

    elif method == 'griddata':
        # NaN 위치만 최근접 보간(불필요한 전점 보간 제거)
        points = np.column_stack((lon2d[valid], lat2d[valid]))
        values = var2d[valid]
        targets = (lon2d[nan_mask], lat2d[nan_mask])
        filled_vals = griddata(points, values, targets, method='nearest')
        var_filled[nan_mask] = filled_vals

    elif method == 'nearestKD':  # 선택지로 추가
        # 지리 최근접(KDTree). 대형 격자에서 빠름.
        # 필요시 지역 투영(예: x=lon*cos(lat0), y=lat)로 왜곡 줄일 수 있음.
        lat0 = np.nanmean(lat2d)
        X = lon2d * np.cos(np.deg2rad(lat0))
        Y = lat2d
        pts_valid = np.column_stack((X[valid], Y[valid]))
        pts_nan   = np.column_stack((X[nan_mask], Y[nan_mask]))
        tree = cKDTree(pts_valid)
        _, idx = tree.query(pts_nan, k=1)
        var_filled[nan_mask] = var2d[valid][idx]

    else:
        raise ValueError(f"Unknown method: {method}")

    # 원 dtype 유지(원하면)
    # return var_filled.astype(var2d.dtype, copy=False)
    return var_filled


def extract_bry(var, direction):
    if direction == 'west':
        return var[..., :, 0]
    elif direction == 'east':
        return var[..., :, -1]
    elif direction == 'south':
        return var[..., 0, :]
    elif direction == 'north':
        return var[..., -1, :]



def get_bottom(var3d, mask, spval=-1e10):
    """
    각 (j, i) column에 대해 bottom-most valid layer index를 반환
    - var3d: shape = (Nz, Ny, Nx)
    - mask:  shape = (Ny, Nx)
    """
    Nz, Ny, Nx = var3d.shape
    bottom = np.zeros((Ny, Nx), dtype=np.int32)

    for j in range(Ny):
        for i in range(Nx):
            if mask[j, i] == 0:
                bottom[j, i] = 0
                continue
            for k in range(Nz - 1, -1, -1):  # 아래부터 위로
                v = var3d[k, j, i]
                if not np.isnan(v) and abs((v - spval)/spval) > 1e-5:
                    bottom[j, i] = k
                    break
            else:
                bottom[j, i] = 0
    return bottom

def get_bottom_vectorized(var3d, mask, spval=-1e10):
    Nz, Ny, Nx = var3d.shape
    valid_mask = (np.abs((var3d - spval) / spval) > 1e-5) & (~np.isnan(var3d))

    bottom = np.zeros((Ny, Nx), dtype=np.int32)
    idx = np.where(mask == 1)

    for j, i in zip(*idx):
        valid_k = np.where(valid_mask[:, j, i])[0]
        if valid_k.size > 0:
            bottom[j, i] = valid_k[-1]
    return bottom

def flood_vertical(varz, mask, spval=-1e10):
    """
    가장 아래 valid 값으로 NaN을 아래쪽에서 채워주는 수직 flood
    - varz: shape = (Nz, Ny, Nx)
    - mask: shape = (Ny, Nx)
    """
    varz = varz.copy()
    Nz, Ny, Nx = varz.shape
    varz[np.isnan(varz)] = spval

    bottom = get_bottom(varz, mask, spval)  
    for j in range(Ny):
        for i in range(Nx):
            if mask[j, i] == 1:
                kbot = bottom[j, i]
                varz[kbot:, j, i] = varz[kbot, j, i]
    return varz
def flood_vertical_vectorized(varz, mask, spval=-1e10):
    varz = varz.copy()
    Nz, Ny, Nx = varz.shape
    varz[np.isnan(varz)] = spval

    bottom = get_bottom_vectorized(varz, mask, spval)

    jj, ii = np.where(mask == 1)
    for j, i in zip(jj, ii):
        k = bottom[j, i]
        if k < Nz:
            varz[k:, j, i] = varz[k, j, i]
    return varz

@jit(nopython=True)
def get_bottom_numba(var3d, mask, spval=-1e10):
    Nz, Ny, Nx = var3d.shape
    bottom = np.zeros((Ny, Nx), dtype=np.int32)

    for j in range(Ny):
        for i in range(Nx):
            if mask[j, i] == 0:
                bottom[j, i] = 0
                continue
            for k in range(Nz - 1, -1, -1):
                v = var3d[k, j, i]
                if v == v and abs((v - spval) / spval) > 1e-5:
                    bottom[j, i] = k
                    break
            else:
                bottom[j, i] = 0
    return bottom


@jit(nopython=True)
def flood_vertical_numba(varz, mask, spval=-1e10):
    Nz, Ny, Nx = varz.shape
    varout = np.copy(varz)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                if varout[k, j, i] != varout[k, j, i]:
                    varout[k, j, i] = spval

    bottom = get_bottom_numba(varout, mask, spval)

    for j in range(Ny):
        for i in range(Nx):
            if mask[j, i] == 1:
                kbot = bottom[j, i]
                val = varout[kbot, j, i]
                for k in range(kbot, Nz):
                    varout[k, j, i] = val

    return varout

def make_all_bry_data_shapes(varnames: list[str], num_steps: int, grd, Ns: int) -> dict[str, dict[str, np.ndarray]]:
    """
    여러 변수에 대해 make_bry_data_shape 호출해서 전체 bry_data 구조 생성
    """
    return {var: make_bry_data_shape(var, num_steps, grd, Ns) for var in varnames}


def make_bry_data_shape(varname: str, num_steps: int, grd, Ns: int):
    """
    벡터/스칼라 변수 구분해서 방향별 bry_data shape 생성기
    병렬 처리 및 staggered grid 고려

    Parameters
    ----------
    varname : str
        변수명 ('temp', 'salt', 'zeta', 'u', 'v', ...)
    num_steps : int
        시간 스텝 수
    grd : ConfigObject
        ROMS 그리드 정보 (lat, lon, mask 포함)
    Ns : int
        수직 레벨 수 (layer_n)

    Returns
    -------
    dict[str, np.ndarray]
        {'west': ..., 'east': ..., 'south': ..., 'north': ...}
    """
    Mp, Lp = grd.lat.shape
    dtype = np.float32

    if varname in ['zeta']:
        shape_map = {
            'west':  (num_steps, Mp),
            'east':  (num_steps, Mp),
            'south': (num_steps, Lp),
            'north': (num_steps, Lp),
        }

    elif varname in ['ubar']:
        shape_map = {
            'west':  (num_steps, Mp),
            'east':  (num_steps, Mp),
            'south': (num_steps, Lp - 1),
            'north': (num_steps, Lp - 1),
        }

    elif varname in ['vbar']:
        shape_map = {
            'west':  (num_steps, Mp - 1),
            'east':  (num_steps, Mp - 1),
            'south': (num_steps, Lp),
            'north': (num_steps, Lp),
        }

    elif varname in ['u']:
        shape_map = {
            'west':  (num_steps, Ns, Mp),
            'east':  (num_steps, Ns, Mp),
            'south': (num_steps, Ns, Lp - 1),
            'north': (num_steps, Ns, Lp - 1),
        }

    elif varname in ['v']:
        shape_map = {
            'west':  (num_steps, Ns, Mp - 1),
            'east':  (num_steps, Ns, Mp - 1),
            'south': (num_steps, Ns, Lp),
            'north': (num_steps, Ns, Lp),
        }

    else:  # 일반 스칼라 변수
        shape_map = {
            'west':  (num_steps, Ns, Mp),
            'east':  (num_steps, Ns, Mp),
            'south': (num_steps, Ns, Lp),
            'north': (num_steps, Ns, Lp),
        }

    return {k: np.zeros(v, dtype=dtype) for k, v in shape_map.items()}















