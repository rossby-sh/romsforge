
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

def vertical_interpolation_parallel2(
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



R1 = vertical_interpolation_parallel2(zr, data, vdepth,
                                     n_jobs=-1,
                                     dedup="mean",
                                     extrap_mode="leading")


R2 = vertical_interpolation_parallel2(zr, data, vdepth,
                                     n_jobs=-1,
                                     dedup="jitter",
                                     extrap_mode="padding",
                                     zsur=0.0,
                                     zbot=-6000.0)  # 너 그리드 최심보다 조금 더 깊게
