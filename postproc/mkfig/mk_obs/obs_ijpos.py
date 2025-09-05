import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from matplotlib.path import Path

def obs_ijpos(grdname, obs_lon, obs_lat, Correction=True, obc_edge=False, strict_griddata=True):
    
    grdname='D:/shjo/ROMS_inputs/NWP4_grd_3_10m_LP.nc'
    nc=Dataset('D:/shjo/ROMS_inputs/obs/roms_obs_phyt_30km_N36_250501_250531_.nc')
    obs_lon = nc['obs_lon'][:]
    obs_lat = nc['obs_lat'][:]
    obc_edge=False; strict_griddata=True
    
    # --- Input shape check (transpose to column vectors) ---
    obs_lon = np.atleast_1d(obs_lon).ravel()
    obs_lat = np.atleast_1d(obs_lat).ravel()

    # --- Read grid file ---
    nc = Dataset(grdname)

    # Variable flags
    got = {
        "spherical": "spherical" in nc.variables,
        "lon_rho": "lon_rho" in nc.variables,
        "lat_rho": "lat_rho" in nc.variables,
        "angle": "angle" in nc.variables,
        "mask_rho": "mask_rho" in nc.variables,
        "coast": "lon_coast" in nc.variables and "lat_coast" in nc.variables,
    }

    spherical = True
    if got["spherical"]:
        spherical_val = nc["spherical"][:]
        if isinstance(spherical_val, np.ndarray) and spherical_val.dtype.char in "SU":
            spherical = spherical_val.tobytes().decode("utf-8").strip().lower() == "t"

    if not got["lon_rho"] or not got["lat_rho"]:
        raise ValueError("Missing lon_rho or lat_rho")

    rlon = nc["lon_rho"][:]
    rlat = nc["lat_rho"][:]
    angle = nc["angle"][:] if got["angle"] else np.zeros_like(rlon)
    rmask = nc["mask_rho"][:] if got["mask_rho"] else np.ones_like(rlon)
    nc.close()

    # --- Define domain polygon ---
    eta, xi = rlon.shape
    Xbox = np.concatenate([
        rlon[0, :],               # 북쪽
        rlon[1:, -1],             # 동쪽
        rlon[-1, -2::-1],         # 남쪽
        rlon[-2::-1, 0]           # 서쪽 ← 여기 수정됨!
    ])
        
    Ybox = np.concatenate([
        rlat[0, :],
        rlat[1:, -1],
        rlat[-1, -2::-1],
        rlat[-2::-1, 0]
    ])

    path = Path(np.vstack((Xbox, Ybox)).T)
    IN = path.contains_points(np.vstack((obs_lon, obs_lat)).T)
    
    bounded = IN.copy() if not obc_edge else IN  # edge 미구현
    
    # from scipy import io
    # bounded_MAT=io.loadmat('D:/shjo/ROMS_inputs/obs/test_bounded.mat')['bounded'].reshape(-1)
    # np.all(bounded==bounded_MAT)
    

    Igrid, Jgrid = np.meshgrid(np.arange(xi), np.arange(eta))
    Igrid = Igrid.astype(float)
    Jgrid = Jgrid.astype(float)
    
    # 마스크 적용
    Igrid[rmask < 1] = np.nan
    Jgrid[rmask < 1] = np.nan
    
    
    if strict_griddata == False :
        # interpolation 위한 준비
        points = np.vstack((rlon[rmask == 1], rlat[rmask == 1])).T
        values_I = Igrid[rmask == 1]
        values_J = Jgrid[rmask == 1]
        
        xi = np.vstack((obs_lon[bounded], obs_lat[bounded])).T
        
        # 최종 위치 보간
        Xgrid = np.full_like(obs_lon, np.nan, dtype=float)
        Ygrid = np.full_like(obs_lat, np.nan, dtype=float)
        Xgrid[bounded] = griddata(points, values_I, xi, method='linear',fill_value=np.nan)
        Ygrid[bounded] = griddata(points, values_J, xi, method='linear',fill_value=np.nan)
    else:

        from scipy.interpolate import griddata
        from scipy.spatial import Delaunay
        
        # --- 보간에 사용할 좌표 평탄화 (Fortran-order로 flatten) ---
        points = np.column_stack((rlon.T.ravel(), rlat.T.ravel()))
        values_I = Igrid.T.ravel()  # 예: I index 값 (또는 X index)
        values_J = Jgrid.T.ravel()  # 예: J index 값 (또는 Y index)
        
        # --- 관측 지점 위치 ---
        xi = np.column_stack((obs_lon[bounded], obs_lat[bounded]))
        
        # --- 삼각형 내부 여부 판단 ---
        tri = Delaunay(points.data)
        mask_inside = tri.find_simplex(xi) >= 0
        
        # --- 결과 초기화 ---
        Xgrid = np.full_like(obs_lon, np.nan, dtype=float)
        Ygrid = np.full_like(obs_lat, np.nan, dtype=float)
        
        # --- 내부 지점만 보간 ---
        interp_I = griddata(points, values_I, xi[mask_inside], method='linear')
        interp_J = griddata(points, values_J, xi[mask_inside], method='linear')
        
        # --- 결과 대입 ---
        idx_valid = np.where(bounded)[0][mask_inside]
        Xgrid[idx_valid] = interp_I
        Ygrid[idx_valid] = interp_J
    

    # --- Final filtering (land mask) ---
    invalid = np.isnan(Xgrid) | np.isnan(Ygrid)
    bounded[invalid] = False

    # --- Optional correction ---
    if Correction:
        # TODO: implement curvilinear correction here if needed
        data_dict = {
            'rlon': rlon,
            'rlat': rlat,
            'angle': angle,
            'obs_lon': obs_lon,
            'obs_lat': obs_lat,
            'bounded': bounded,
            'Xgrid': Xgrid,
            'Ygrid': Ygrid
        }
        
        np.save('D:/shjo/ROMS_inputs/obs/correction_inputs.npy', data_dict)

        Xgrid1, Ygrid1 = correct_fractional_coordinates1(rlon.T, rlat.T, angle.T, obs_lon, obs_lat, bounded, Xgrid, Ygrid)
        Xgrid2, Ygrid2 = correct_fractional_coordinates2(rlon, rlat, angle, obs_lon, obs_lat, bounded, Xgrid, Ygrid)

        pass

    return Xgrid, Ygrid


def correct_fractional_coordinates1(rlon, rlat, angle, obs_lon, obs_lat, bounded, X, Y, debug=True):
    """
    Apply curvilinear coordinate correction to fractional X, Y locations.
    """
    import numpy as np
    from netCDF4 import Dataset
    from scipy.interpolate import griddata
    from matplotlib.path import Path
    from scipy  import io

    
    Nobs = len(obs_lon)               # number of observations
    
    block_length = 100000            # size of block arrays
    
    Eradius = 6371315.0              # Earth radius (meters)
    deg2rad = np.pi / 180.0          # degrees to radians factor
    
    Lr, Mr = rlon.shape              # Number of rho-points
    
    Lm = Lr - 1
    Mm = Mr - 1
    
    # --------------------------------------------------------------------------
    #  Set (I,J) coordinates of the grid cell containing the observation
    #  need to add 1 because zero lower bound in ROMS "rlon"
    # --------------------------------------------------------------------------
    # I, J를 -1로 초기화 (shape 유지)
    I = np.full_like(X, fill_value=-1, dtype=int)
    J = np.full_like(Y, fill_value=-1, dtype=int)
    
    # NaN이 아닌 위치만 계산
    valid_X = ~np.isnan(X)
    valid_Y = ~np.isnan(Y)
    
    I[valid_X] = np.floor(X[valid_X]).astype(int)
    J[valid_Y] = np.floor(Y[valid_Y]).astype(int)
    
    # rlon shape
    Lr, Mr = rlon.shape
    Lm = Lr - 1
    Mm = Mr - 1
    
    # 조건 만족하는 위치에 대해서만 +1
    ind_I = np.where((I >= 0) & (I < Lm))[0]
    I[ind_I] += 1
    
    ind_J = np.where((J >= 0) & (J < Mm))[0]
    J[ind_J] += 1
    
    if debug:
        print()
        print(f"  Xmin = {np.nanmin(X):6.2f}  Xmax = {np.nanmax(X):6.2f}  "
              f"Imin = {np.min(I):3d}  Imax = {np.max(I):3d}  Lr = {Lr:3d}")
        
        print(f"  Ymin = {np.nanmin(Y):6.2f}  Ymax = {np.nanmax(Y):6.2f}  "
              f"Jmin = {np.min(J):3d}  Jmax = {np.max(J):3d}  Mr = {Mr:3d}")
        print()
    
    # --------------------------------------------------------------------------
    #  It is possible that we are processing a large number of observations.
    #  Therefore, the observation vector is processed by blocks to reduce
    #  the memory requirements.
    # --------------------------------------------------------------------------
    
    N = int(np.ceil(Nobs / block_length))
    n1 = 0
    n2 = 0
    
    
    while n2 < Nobs:
        n1 = n2 + 1
        n2 = n1 + block_length
        if n2 > Nobs:
            n2 = Nobs
    
        # 블록 인덱스: [n1, n2) ← 파이썬은 end 포함 안 함
        iobs = np.arange(n1 - 1, n2)
    
        # bounded = boolean array of length Nobs
        ind = np.where(~bounded.astype('bool')[n1-1:n2])[0]
        if ind.size > 0:
            iobs = np.delete(iobs, ind)
    
        # I, J는 전체 길이 Nobs에 대해 미리 계산된 int array
        # rlon.shape = (Lr, Mr) = (ny, nx)
        shape = rlon.shape
        ny, nx = shape
    
    
        i_j   = np.ravel_multi_index((I[iobs]-1,     J[iobs]-1),     dims=shape, order='F')
        ip1_j = np.ravel_multi_index((I[iobs]-1 + 1, J[iobs]-1),     dims=shape, order='F')
        i_jp1 = np.ravel_multi_index((I[iobs]-1,     J[iobs]-1 + 1), dims=shape, order='F')
    
        if debug:
            print(f"  Processing observation vector, n1:n2 = {n1:07d} - {n2:07d}  size = {len(iobs)}")
        
        
        # 단위 변환 계수
        yfac = Eradius * deg2rad
        xfac = yfac * np.cos(obs_lat[iobs] * deg2rad)
        
        rlon_flat = rlon.flatten(order='F')
        rlat_flat = rlat.flatten(order='F')
        
        xpp = (obs_lon[iobs] - rlon_flat[i_j]) * xfac
        ypp = (obs_lat[iobs] - rlat_flat[i_j]) * yfac
        
    
        # Law of Cosines로 셀 비틀림 각도 계산
        diag2 = (rlon_flat[ip1_j] - rlon_flat[i_jp1])**2 + \
                (rlat_flat[ip1_j] - rlat_flat[i_jp1])**2
        
        aa2 = (rlon_flat[i_j] - rlon_flat[ip1_j])**2 + \
              (rlat_flat[i_j] - rlat_flat[ip1_j])**2
        
        bb2 = (rlon_flat[i_j] - rlon_flat[i_jp1])**2 + \
              (rlat_flat[i_j] - rlat_flat[i_jp1])**2
        
        phi = np.arcsin((diag2 - aa2 - bb2) / (2.0 * np.sqrt(aa2 * bb2)))
        
        
        angle_flat = angle.flatten(order='F')
        
        # --- 1. curvilinear 회전 적용 ---
        dx = xpp * np.cos(angle_flat[i_j]) + ypp * np.sin(angle_flat[i_j])
        dy = ypp * np.cos(angle_flat[i_j]) - xpp * np.sin(angle_flat[i_j])
        
        # --- 2. 평행사변형 보정 ---
        dx = dx + dy * np.tan(phi)
        dy = dy / np.cos(phi)
        
        # --- 3. 셀 길이로 정규화해서 fractional dx로 변환 ---
        dx = dx / np.sqrt(aa2) / xfac
        
        dx = np.clip(dx, 0, 1)
        
        # dy는 스케일 변환이 추가됨
        dy /= np.sqrt(bb2) * yfac
        dy = np.clip(dy, 0, 1)
        
        # 최종 fractional 좌표 계산
        X[iobs] = np.floor(X[iobs]) + dx
        Y[iobs] = np.floor(Y[iobs]) + dy
    
        return X, Y



def correct_fractional_coordinates2(rlon, rlat, angle, obs_lon, obs_lat, bounded, X, Y, debug=True):
    import numpy as np

    debugging = debug

    Nobs = len(obs_lon)
    block_length = 100000
    Eradius = 6371315.0
    deg2rad = np.pi / 180.0

    ny, nx = rlon.shape
    Lm, Mm = nx - 1, ny - 1

    I = np.full_like(X, fill_value=-1, dtype=int)
    J = np.full_like(Y, fill_value=-1, dtype=int)

    valid_X = ~np.isnan(X)
    valid_Y = ~np.isnan(Y)
    I[valid_X] = np.floor(X[valid_X]).astype(int)
    J[valid_Y] = np.floor(Y[valid_Y]).astype(int)

    I[(I >= 0) & (I < Lm)] += 1
    J[(J >= 0) & (J < Mm)] += 1

    if debugging:
        print()
        print(f"  Xmin = {np.nanmin(X):6.2f}  Xmax = {np.nanmax(X):6.2f}  Imin = {np.min(I):3d}  Imax = {np.max(I):3d}  nx = {nx:3d}")
        print(f"  Ymin = {np.nanmin(Y):6.2f}  Ymax = {np.nanmax(Y):6.2f}  Jmin = {np.min(J):3d}  Jmax = {np.max(J):3d}  ny = {ny:3d}")
        print()

    N = int(np.ceil(Nobs / block_length))
    n1, n2 = 0, 0

    while n2 < Nobs:
        n1 = n2
        n2 = min(n1 + block_length, Nobs)
        iobs = np.arange(n1, n2)

        ind = np.where(~bounded[iobs])[0]
        if ind.size > 0:
            iobs = np.delete(iobs, ind)

        shape = (ny, nx)
        i_j   = np.ravel_multi_index((J[iobs]-1,     I[iobs]-1), shape, order='C')
        ip1_j = np.ravel_multi_index((J[iobs]-1,     I[iobs]),   shape, order='C')
        i_jp1 = np.ravel_multi_index((J[iobs],       I[iobs]-1), shape, order='C')

        yfac = Eradius * deg2rad
        xfac = yfac * np.cos(obs_lat[iobs] * deg2rad)

        rlon_flat = rlon.flatten(order='C')
        rlat_flat = rlat.flatten(order='C')

        xpp = (obs_lon[iobs] - rlon_flat[i_j]) * xfac
        ypp = (obs_lat[iobs] - rlat_flat[i_j]) * yfac

        diag2 = (rlon_flat[ip1_j] - rlon_flat[i_jp1])**2 + (rlat_flat[ip1_j] - rlat_flat[i_jp1])**2
        aa2 = (rlon_flat[i_j] - rlon_flat[ip1_j])**2 + (rlat_flat[i_j] - rlat_flat[ip1_j])**2
        bb2 = (rlon_flat[i_j] - rlon_flat[i_jp1])**2 + (rlat_flat[i_j] - rlat_flat[i_jp1])**2

        phi = np.arcsin((diag2 - aa2 - bb2) / (2.0 * np.sqrt(aa2 * bb2)))

        angle_flat = angle.flatten(order='C')
        dx = xpp * np.cos(angle_flat[i_j]) + ypp * np.sin(angle_flat[i_j])
        dy = ypp * np.cos(angle_flat[i_j]) - xpp * np.sin(angle_flat[i_j])

        dx += dy * np.tan(phi)
        dy /= np.cos(phi)

        dx = dx / np.sqrt(aa2) / xfac
        dy = dy / (np.sqrt(bb2) * yfac)

        dx = np.clip(dx, 0, 1)
        dy = np.clip(dy, 0, 1)

        X[iobs] = np.floor(X[iobs]) + dx
        Y[iobs] = np.floor(Y[iobs]) + dy

    return X, Y




















