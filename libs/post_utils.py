import numpy as np
from joblib import Parallel, delayed
from scipy.interpolate import interp1d


# ---------- helpers ----------

def _fill_surface_bottom_only(col: np.ndarray) -> np.ndarray:
    """leading/trailing NaN만 최근접 유효값으로 채움(내부 NaN은 그대로)."""
    out = col.copy()
    finite = np.isfinite(out)
    if not finite.any():
        return out
    # 표층(앞쪽)
    first = np.argmax(finite)
    if first > 0:
        out[:first] = out[first]
    # 저층(뒤쪽)
    last = len(out) - 1 - np.argmax(finite[::-1])
    if last < len(out) - 1:
        out[last + 1:] = out[last]
    return out


def _enforce_monotonic(z: np.ndarray,
                       v: np.ndarray,
                       dedup: str = "jitter",
                       eps: float = 1e-8) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """
    NaN 제거 → z 오름차순 정렬(얕은→깊은) → 중복 처리.
    dedup: 'jitter' | 'mean' | 'none'
    """
    m = np.isfinite(z) & np.isfinite(v)
    if m.sum() < 2:
        return None, None
    z = z[m].astype(float)
    v = v[m].astype(float)

    # 오름차순 정렬(얕은→깊은). (깊이가 음수라도 '값 오름차순'으로 맞춤)
    order = np.argsort(z)
    z = z[order]
    v = v[order]

    if dedup == "jitter":
        # 거의-중복/동일 z에 미세 증가분 부여
        dz = np.diff(z)
        k = np.where(dz <= 0)[0]
        if k.size:
            # 스케일 적응 eps
            scale = max(abs(z[-1] - z[0]), 1.0)
            step = eps * scale
            for i in k:
                if z[i + 1] <= z[i]:
                    z[i + 1] = z[i] + step
    elif dedup == "mean":
        zu, inv = np.unique(z, return_inverse=True)
        if len(zu) != len(z):
            acc = np.zeros_like(zu, dtype=float)
            cnt = np.zeros_like(zu, dtype=int)
            for t, val in zip(inv, v):
                acc[t] += val; cnt[t] += 1
            z = zu
            v = acc / np.maximum(cnt, 1)
    # 'none'이면 그대로

    if z.size < 2:
        return None, None
    return z, v


def _sort_stdz(stdvdepth: np.ndarray) -> tuple[np.ndarray, np.ndarray | None]:
    """
    stdvdepth가 증가/감소 어떤 단조든 상관없이,
    내부 보간은 '오름차순'으로 진행하고, 나중에 원순서 복원.
    """
    stdz = np.asarray(stdvdepth, dtype=float)
    if stdz.size == 0:
        return stdz, None
    # 이미 오름차순이면 그대로
    if np.all(np.diff(stdz) > 0):
        return stdz, None
    order = np.argsort(stdz)
    return stdz[order], np.argsort(order)


# ---------- main ----------

def interpolate_s_to_zlevels(
    depth3d: np.ndarray,
    roms_var: np.ndarray,
    stdvdepth: np.ndarray,
    *,
    n_jobs: int = -1,
    dedup: str = "jitter",            # 'jitter' | 'mean' | 'none'
    extrap_mode: str = "leading",     # 'leading' | 'padding'
    zsur: float | None = None,        # extrap_mode='padding'일 때 필요(보통 0.0)
    zbot: float | None = None         # extrap_mode='padding'일 때 필요(충분히 깊은 음수)
) -> np.ndarray:
    """
    s-레벨 → z-레벨 수직 보간.
    - extrap_mode='leading' : roms2z 스타일 (보간은 경계 NaN, 마지막에 표층/저층만 메움)
    - extrap_mode='padding' : interp_field 스타일 (보간 전 경계 패딩으로 최근접 외삽 강제)
    """
    Ns, Ny, Nx = depth3d.shape
    Nz = len(stdvdepth)

    # stdvdepth을 내부 오름차순으로 정렬, 나중에 복원
    stdz_sorted, restore = _sort_stdz(stdvdepth)

    if extrap_mode not in ("leading", "padding"):
        raise ValueError("extrap_mode는 'leading' 또는 'padding'만 지원")

    if extrap_mode == "padding" and (zsur is None or zbot is None):
        raise ValueError("extrap_mode='padding'일 때 zsur, zbot 반드시 지정")

    def interp_column(ix: int) -> np.ndarray:
        out = np.full((Nz, Ny), np.nan, dtype=float)
        for jy in range(Ny):
            z_col = depth3d[:, jy, ix]
            v_col = roms_var[:, jy, ix]

            z, v = _enforce_monotonic(z_col, v_col, dedup=dedup)
            if z is None:
                continue

            if extrap_mode == "padding":
                # interp_field (rutergus) 스타일: 경계 패딩 후 보간(항상 최근접 외삽)
                z_pad = np.concatenate(([zbot], z, [zsur]))
                v_pad = np.concatenate(([v[0]], v, [v[-1]]))
                f = interp1d(z_pad, v_pad, kind="linear",
                             bounds_error=False, fill_value="extrapolate",
                             assume_sorted=True)
                col = f(stdz_sorted)
            else:
                # roms2z (pyroms) 스타일: 경계는 NaN → 나중에 표층/저층만 채움
                f = interp1d(z, v, kind="linear",
                             bounds_error=False, fill_value=(np.nan, np.nan),
                             assume_sorted=True)
                col = f(stdz_sorted)
                col = _fill_surface_bottom_only(col)

            # 원 stdvdepth 순서 복원
            if restore is not None:
                tmp = np.empty_like(col)
                tmp[restore] = col
                col = tmp

            out[:, jy] = col
        return out

    cols = Parallel(n_jobs=n_jobs, prefer="threads")(delayed(interp_column)(i) for i in range(Nx))
    std_area = np.stack(cols, axis=-1)  # (Nz, Ny, Nx)
    return std_area




import numpy as np
from joblib import Parallel, delayed
from scipy.interpolate import interp1d

def vertical_interp_to_ZR(
    ZD: np.ndarray,          # (Kd, Ny, Nx) : donor depth on receiver XY (수평 리맵된 z_donor)
    VD: np.ndarray,          # (Kd, Ny, Nx) : donor var on receiver XY (수평 리맵된 var)
    ZR: np.ndarray,          # (Kr, Ny, Nx) : receiver depth (타깃 z_rho 또는 z_u/z_v)
    *,
    n_jobs: int = -1,
    dedup: str = "jitter",   # 'jitter' | 'mean' | 'none'  (거의-중복 z 처리)
    extrap_mode: str = "padding",   # 'padding' | 'leading'  (경계 외삽 방식)
    zsur: float | None = None,  # padding 모드: None이면 컬럼별 max(ZD) 사용
    zbot: float | None = None,  # padding 모드: None이면 컬럼별 min(ZD) 사용
    kind: str = "linear",    # 'linear'|'nearest'|'cubic'...
):
    def _to_float_nan(a):
        return a.filled(np.nan).astype(float, copy=False) if np.ma.isMaskedArray(a) else np.asarray(a, dtype=float)

    def _fill_surface_bottom_only(col):
        out = col.copy()
        finite = np.isfinite(out)
        if not finite.any():
            return out
        first = np.argmax(finite)
        if first > 0: out[:first] = out[first]
        last = len(out) - 1 - np.argmax(finite[::-1])
        if last < len(out) - 1: out[last+1:] = out[last]
        return out

    def _enforce_monotonic(z, v, dedup="jitter", eps=1e-8):
        m = np.isfinite(z) & np.isfinite(v)
        if m.sum() < 2: return None, None
        z = z[m].astype(float, copy=False); v = v[m].astype(float, copy=False)
        o = np.argsort(z); z = z[o]; v = v[o]
        if dedup == "jitter":
            dz = np.diff(z); k = np.where(dz <= 0)[0]
            if k.size:
                scale = max(abs(z[-1]-z[0]), 1.0); step = eps*scale
                for i in k:
                    if z[i+1] <= z[i]: z[i+1] = z[i] + step
        elif dedup == "mean":
            zu, inv = np.unique(z, return_inverse=True)
            if len(zu) != len(z):
                acc = np.zeros_like(zu, float); cnt = np.zeros_like(zu, int)
                for t, val in zip(inv, v): acc[t]+=val; cnt[t]+=1
                z = zu; v = acc/np.maximum(cnt,1)
        if z.size < 2: return None, None
        return z, v

    ZD = _to_float_nan(ZD); VD = _to_float_nan(VD); ZR = _to_float_nan(ZR)
    Kd, Ny, Nx = ZD.shape
    Kr = ZR.shape[0]
    assert VD.shape == (Kd, Ny, Nx) and ZR.shape[1:] == (Ny, Nx)

    def interp_col(ix):
        out = np.full((Kr, Ny), np.nan, float)
        for jy in range(Ny):
            zd = ZD[:, jy, ix]
            vd = VD[:, jy, ix]
            zr = ZR[:, jy, ix]

            z, v = _enforce_monotonic(zd, vd, dedup=dedup)
            if z is None:
                continue

            if extrap_mode == "padding":
                zbot_here = np.nanmin(z) if zbot is None else float(zbot)
                zsur_here = np.nanmax(z) if zsur is None else float(zsur)
                z_pad = np.concatenate(([zbot_here], z, [zsur_here]))
                v_pad = np.concatenate(([v[0]],      v, [v[-1]     ]))
                f = interp1d(z_pad, v_pad, kind=kind, bounds_error=False,
                             fill_value="extrapolate", assume_sorted=True)
                col = f(zr)
            else:
                f = interp1d(z, v, kind=kind, bounds_error=False,
                             fill_value=(np.nan, np.nan), assume_sorted=True)
                col = f(zr)
                col = _fill_surface_bottom_only(col)

            out[:, jy] = col
        return out

    cols = Parallel(n_jobs=n_jobs, prefer="threads")(delayed(interp_col)(i) for i in range(Nx))
    return np.stack(cols, axis=-1)  # (Kr, Ny, Nx)




import numpy as np
from numba import njit, prange

@njit
def _fill_surface_bottom_only(col):
    n = col.shape[0]
    # first finite
    first = -1
    for i in range(n):
        if not np.isnan(col[i]):
            first = i
            break
    if first == -1:
        return
    for i in range(first):
        col[i] = col[first]
    # last finite
    last = -1
    for i in range(n - 1, -1, -1):
        if not np.isnan(col[i]):
            last = i
            break
    for i in range(last + 1, n):
        col[i] = col[last]

@njit
def _make_strictly_increasing(z):
    n = z.shape[0]
    if n <= 1:
        return
    span = abs(z[n - 1] - z[0])
    if span < 1.0:
        span = 1.0
    step = 1e-8 * span
    for i in range(n - 1):
        if z[i + 1] <= z[i]:
            z[i + 1] = z[i] + step

@njit
def _linear_interp_sorted(x, y, xq, out, extrapolate):
    """
    x: strictly increasing
    y: same length
    xq: query points
    out: preallocated
    extrapolate: True → linear extrap at both ends, False → NaN at both ends
    """
    n = x.shape[0]
    m = xq.shape[0]
    for j in range(m):
        q = xq[j]
        # left
        if q <= x[0]:
            if extrapolate:
                if n >= 2:
                    dx = x[1] - x[0]
                    if dx == 0.0:
                        out[j] = y[0]
                    else:
                        t = (q - x[0]) / dx
                        out[j] = y[0] + t * (y[1] - y[0])
                else:
                    out[j] = y[0]
            else:
                out[j] = np.nan
            continue
        # right
        if q >= x[n - 1]:
            if extrapolate:
                if n >= 2:
                    dx = x[n - 1] - x[n - 2]
                    if dx == 0.0:
                        out[j] = y[n - 1]
                    else:
                        t = (q - x[n - 2]) / dx
                        out[j] = y[n - 2] + t * (y[n - 1] - y[n - 2])
                else:
                    out[j] = y[n - 1]
            else:
                out[j] = np.nan
            continue
        # binary search (left)
        lo = 0
        hi = n
        while lo < hi:
            mid = (lo + hi) // 2
            if x[mid] < q:
                lo = mid + 1
            else:
                hi = mid
        i = lo - 1
        x0 = x[i]; x1 = x[i + 1]
        y0 = y[i]; y1 = y[i + 1]
        dx = x1 - x0
        if dx == 0.0:
            out[j] = 0.5 * (y0 + y1)
        else:
            t = (q - x0) / dx
            out[j] = y0 + t * (y1 - y0)

@njit(parallel=True)
def vertical_interp_to_ZR_numba(ZD, VD, ZR, mode_padding=True, clamp=True):
    """
    ZD: (Kd, Ny, Nx)  donor depth on receiver XY (already horizontally remapped)
    VD: (Kd, Ny, Nx)  donor variable on receiver XY (already horizontally remapped)
    ZR: (Kr, Ny, Nx)  receiver target depths (z_rho or z_u/z_v of target grid)

    mode_padding=True  → MATLAB interp_field 스타일 (양끝 선형 외삽)
    mode_padding=False → leading 스타일 (경계 NaN, 이후 표층/저층만 채움)
    clamp=True         → 결과를 각 칼럼의 [min(VD_col), max(VD_col)]로 제한

    return: (Kr, Ny, Nx) array
    """
    Kd, Ny, Nx = ZD.shape
    Kr = ZR.shape[0]
    out = np.empty((Kr, Ny, Nx), dtype=np.float64)

    for j in prange(Ny):
        for i in range(Nx):
            # collect finite pairs
            cnt = 0
            for k in range(Kd):
                z = ZD[k, j, i]
                v = VD[k, j, i]
                if not (np.isnan(z) or np.isnan(v)):
                    cnt += 1
            if cnt < 2:
                for kk in range(Kr):
                    out[kk, j, i] = np.nan
                continue

            zc = np.empty(cnt, dtype=np.float64)
            vc = np.empty(cnt, dtype=np.float64)
            p = 0
            vmin = 1e300
            vmax = -1e300
            for k in range(Kd):
                z = ZD[k, j, i]
                v = VD[k, j, i]
                if not (np.isnan(z) or np.isnan(v)):
                    zc[p] = z
                    vc[p] = v
                    if v < vmin: vmin = v
                    if v > vmax: vmax = v
                    p += 1

            # sort by depth (ascending)
            for a in range(1, cnt):
                zt = zc[a]; vt = vc[a]; b = a - 1
                while b >= 0 and zc[b] > zt:
                    zc[b + 1] = zc[b]; vc[b + 1] = vc[b]; b -= 1
                zc[b + 1] = zt; vc[b + 1] = vt
            _make_strictly_increasing(zc)

            # interpolate
            col = np.empty(Kr, dtype=np.float64)
            _linear_interp_sorted(zc, vc, ZR[:, j, i], col, extrapolate=mode_padding)

            # leading: fill only surface/bottom
            if not mode_padding:
                _fill_surface_bottom_only(col)

            # clamp to donor column range
            if clamp:
                for kk in range(Kr):
                    vq = col[kk]
                    if not np.isnan(vq):
                        if vq < vmin: vq = vmin
                        elif vq > vmax: vq = vmax
                        col[kk] = vq

            for kk in range(Kr):
                out[kk, j, i] = col[kk]

    return out

def vertical_interp_to_ZR2(ZD, VD, ZR, *, mode="padding", clamp=True):
    """
    깔끔한 래퍼: 마스크/dtype 정리 후 numba 함수 호출
    mode: 'padding' | 'leading'
    """
    def _to_float(a):
        if np.ma.isMaskedArray(a):
            return a.filled(np.nan).astype(np.float64, copy=False)
        return np.asarray(a, dtype=np.float64)

    ZD_ = _to_float(ZD)
    VD_ = _to_float(VD)
    ZR_ = _to_float(ZR)

    assert ZD_.shape == VD_.shape, "ZD and VD must have the same shape (Kd, Ny, Nx)"
    assert ZR_.shape[1:] == ZD_.shape[1:], "ZR (Kr, Ny, Nx) must match Ny,Nx of ZD/VD"

    mode_padding = (mode == "padding")
    return vertical_interp_to_ZR_numba(ZD_, VD_, ZR_, mode_padding=mode_padding, clamp=clamp)





def _pad_nearest(a, axis):
    sl_lo = [slice(None)]*a.ndim
    sl_hi = [slice(None)]*a.ndim
    sl_lo[axis] = 0
    sl_hi[axis] = -1
    lo = a[tuple(sl_lo)]
    hi = a[tuple(sl_hi)]
    return np.concatenate([np.expand_dims(lo, axis), a, np.expand_dims(hi, axis)], axis=axis)

# def _nan_pair_mean(a, b):
#     w1 = np.isfinite(a).astype(np.float32)
#     w2 = np.isfinite(b).astype(np.float32)
#     num = np.where(np.isfinite(a), a, 0.0) * w1 + np.where(np.isfinite(b), b, 0.0) * w2
#     den = w1 + w2
#     out = np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den>0)
#     return out


def _nan_pair_mean(a, b, *, keep_nan_if_any=True):
    """
    keep_nan_if_any=True  -> 두 값 중 하나라도 NaN이면 결과 NaN (엄격)
    keep_nan_if_any=False -> 기존처럼 유효한 것들만 평균
    """
    a_fin = np.isfinite(a)
    b_fin = np.isfinite(b)

    if keep_nan_if_any:
        both = a_fin & b_fin
        out = np.full(a.shape, np.nan, dtype=np.result_type(a, b, float))
        out[both] = 0.5*(a[both] + b[both])
        return out
    else:
        w1 = a_fin.astype(np.float32)
        w2 = b_fin.astype(np.float32)
        num = np.where(a_fin, a, 0.0) + np.where(b_fin, b, 0.0)
        den = w1 + w2
        out = np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den>0)
        return out

def u2rho_rutgers_safenan(u):
    """
    u: (..., eta, xi-1) -> rho: (..., eta, xi)
    Rutgers(가장자리 복제) + 인접 2점 NaN-무시 평균
    """
    if u.ndim == 2:
        u = u[np.newaxis, ...]
    up = _pad_nearest(u, axis=-1)                # (..., eta, xi+1)
    rho = _nan_pair_mean(up[..., :, :-1], up[..., :, 1:])  # 1칸 쉬프트 평균
    return np.squeeze(rho)

def v2rho_rutgers_safenan(v):
    """
    v: (..., eta-1, xi) -> rho: (..., eta, xi)
    rutgers(가장자리 복제) + 인접 2점 nan-무시 평균
    """
    if v.ndim == 2:
        v = v[np.newaxis, ...]
    vp = _pad_nearest(v, axis=-2)                # (..., eta+1, xi)
    rho = _nan_pair_mean(vp[..., :-1, :], vp[..., 1:, :])  # 1칸 쉬프트 평균
    return np.squeeze(rho)


def uv2rho_rutgers_safenan(var,uv="u"):
    """
    u: (..., eta, xi-1) -> rho: (..., eta, xi)
    Rutgers(가장자리 복제) + 인접 2점 NaN-무시 평균
    """
    if uv=="u":
        if var.ndim == 2:
            var = var[np.newaxis, ...]
        up = _pad_nearest(var, axis=-1)                # (..., eta, xi+1)
        rho = _nan_pair_mean(up[..., :, :-1], up[..., :, 1:])  # 1칸 쉬프트 평균
        return np.squeeze(rho)
    elif uv=="v":
        if var.ndim == 2:
            var = var[np.newaxis, ...]
        vp = _pad_nearest(var, axis=-2)                # (..., eta+1, xi)
        rho = _nan_pair_mean(vp[..., :-1, :], vp[..., 1:, :])  # 1칸 쉬프트 평균
        return np.squeeze(rho)



# ---------- 2) 라플라시안 스무딩 인페인팅 (NaN 영역 메움) ----------
def inpaint_nan_smooth(arr, mask=None, max_iter=300, tol=1e-5):
    """
    2D 또는 3D 배열의 NaN을 4-이웃 평균으로 반복 메움(라플라시안 완화).
    - mask: rho 마스크(1=바다, 0=육지). 주면 육지(0)는 절대 채우지 않음.
    - 3D는 배치 축(첫 축)마다 독립 처리.
    """
    def _inpaint2d(a2d, m2d):
        a = a2d.copy()
        nan0 = ~np.isfinite(a)
        if m2d is not None:
            # 육지는 고정 경계로 취급 (채우지 않음)
            nan_fill = nan0 & (m2d.astype(bool))   # 바다에서만 채움
            fixed = (~nan_fill)
        else:
            nan_fill = nan0
            fixed = ~nan_fill

        if not np.any(nan_fill):
            return a

        # 초기값: 주변 유효값의 평균으로 빠르게 시드
        # (없으면 0으로 시작)
        north = np.roll(a, 1, axis=0)
        south = np.roll(a, -1, axis=0)
        west  = np.roll(a, 1, axis=1)
        east  = np.roll(a, -1, axis=1)
        w = np.isfinite(north).astype(float) + np.isfinite(south).astype(float) + \
            np.isfinite(west).astype(float)  + np.isfinite(east).astype(float)
        seed = (
            np.where(np.isfinite(north), north, 0.0) +
            np.where(np.isfinite(south), south, 0.0) +
            np.where(np.isfinite(west),  west,  0.0) +
            np.where(np.isfinite(east),  east,  0.0)
        )
        with np.errstate(invalid='ignore'):
            a[nan_fill] = np.where(w[nan_fill] > 0, seed[nan_fill] / w[nan_fill], 0.0)

        # Gauss-Seidel 반복
        for _ in range(max_iter):
            old = a.copy()

            # 4-이웃 평균
            north = np.roll(a, 1, axis=0)
            south = np.roll(a, -1, axis=0)
            west  = np.roll(a, 1, axis=1)
            east  = np.roll(a, -1, axis=1)

            wN = np.isfinite(north).astype(float)
            wS = np.isfinite(south).astype(float)
            wW = np.isfinite(west ).astype(float)
            wE = np.isfinite(east ).astype(float)
            wsum = wN + wS + wW + wE
            newv = (
                np.where(np.isfinite(north), north, 0.0) +
                np.where(np.isfinite(south), south, 0.0) +
                np.where(np.isfinite(west ), west,  0.0) +
                np.where(np.isfinite(east ), east,  0.0)
            )
            with np.errstate(invalid='ignore'):
                a[nan_fill] = np.where(wsum[nan_fill] > 0, newv[nan_fill]/wsum[nan_fill], a[nan_fill])

            # 고정점(비-채움 영역)은 유지
            a[~np.isfinite(a) & fixed] = old[~np.isfinite(a) & fixed]

            # 수렴 체크
            delta = np.nanmax(np.abs(a - old))
            if not np.isfinite(delta) or delta < tol:
                break
        return a

    arr = np.asarray(arr)
    if arr.ndim == 2:
        return _inpaint2d(arr, mask)
    elif arr.ndim == 3:
        out = []
        for k in range(arr.shape[0]):
            out.append(_inpaint2d(arr[k], mask))
        return np.stack(out, axis=0)
    else:
        raise ValueError("arr must be 2D or 3D")


