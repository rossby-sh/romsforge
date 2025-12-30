
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
avg_3d_monthly_avg_ol_nc4.py
- 오버랩(4일) 스킵 + 전층(3D) 유지 + 월평균 + 15일 정규화
- xarray/dask 없이 netCDF4.Dataset으로 순차 스트리밍
- 파일명 끝 '<cycle:4d>_<chunk:4d>.nc'만 파싱 (접두어 변화 무관)
- 첫 사이클: 0001~0008, 이후 사이클: 0005~0008만 사용
"""

# --------- NFS/HDF5 안전 패치 ---------
import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# ---------------- 임포트 ----------------
import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from netCDF4 import Dataset as NC, num2date, date2num

# --- libs/utils.py 로드 ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl  # noqa: E402


# ---------------- [01] 설정 로드 ----------------
cfg = tl.parse_config("./config_proc_nft.yaml")
_ = tl.load_roms_grid(getattr(cfg, "grdname", None))  # 필요 시 확인용

case = getattr(cfg, "case", "case")
data_dir = Path(getattr(cfg, "data_dir", ".")).expanduser().resolve()
base_dir = Path(getattr(cfg, "base_dir", ".")).expanduser().resolve()
outdir = base_dir / str(case)
outdir.mkdir(parents=True, exist_ok=True)

data_prefix_header = str(getattr(cfg, "data_prefix_header", ""))  # 비워도 허용
VARS = list(getattr(cfg, "vars", ["zeta", "u", "v", "temp", "salt"]))
chunks_per_cycle = int(getattr(cfg, "chunks_per_cycle", 8))
overlap = int(getattr(cfg, "overlap", 4))
validate = bool(getattr(cfg, "validate", True))

vn = getattr(cfg, "var_names", None) or {}
time_candidates = list(getattr(vn, "time_candidates", ["ocean_time", "time", "t"]))
ref_time_name   = getattr(vn, "ref_time", None)


# ---------------- [02] 입력 파일 스캔 ----------------
cands = []
for p in data_dir.iterdir():
    if not p.is_file():
        continue
    if p.suffix.lower() not in {".nc", ".nc4", ".cdf"}:
        continue
    if data_prefix_header and not p.name.startswith(data_prefix_header):
        continue
    cands.append(str(p))
files_all = sorted(cands)
if not files_all:
    raise FileNotFoundError(f"입력 파일 없음: {data_dir}, prefix={data_prefix_header!r}")

print(f"[info] files={len(files_all)}개 (예: {files_all[:3]})")
print(f"[info] 평균 변수: {VARS}")


# ---------------- [03] 오버랩 제거 ----------------
NAME_RE = re.compile(r"(?P<cycle>\d{4})_(?P<chunk>\d{4})\.nc$")

def select_nonoverlap(all_files, chunks_per_cycle=8, overlap=4):
    groups = {}  # cycle -> [(chunk, path)]
    for f in all_files:
        m = NAME_RE.search(Path(f).name)
        if not m:
            continue
        cyc = int(m.group("cycle")); chk = int(m.group("chunk"))
        groups.setdefault(cyc, []).append((chk, f))
    if not groups:
        return all_files
    cycles = sorted(groups)
    first_cycle = cycles[0]
    chosen = []
    for cyc in cycles:
        s = 1 if cyc == first_cycle else overlap + 1
        e = chunks_per_cycle
        for chk, f in sorted(groups[cyc]):
            if s <= chk <= e:
                chosen.append(f)
    return chosen

files = select_nonoverlap(files_all, chunks_per_cycle=chunks_per_cycle, overlap=overlap)
print(f"[info] 비중복 선택 후 files={len(files)}개")

# 디버그: 사이클별 선택 개수 요약
cnt = defaultdict(int)
for f in files:
    m = NAME_RE.search(Path(f).name)
    if m:
        cnt[int(m.group("cycle"))] += 1
if cnt:
    sample = sorted(cnt.items())[:10]
    mean_sel = np.mean([v for _, v in cnt.items()])
    print(f"[check] 사이클별 선택 개수(앞5): {sample}  (평균≈{mean_sel:.2f})")


# ---------------- [04] 유효 파일 검증(옵션) ----------------
def filter_valid_netcdf(paths):
    if not validate:
        return paths
    ok, bad = [], []
    for f in paths:
        try:
            with NC(f, "r"):
                pass
            ok.append(f)
        except Exception as e:
            bad.append((f, str(e)))
    if bad:
        print(f"[warn] 손상/비호환 파일 {len(bad)}개 제외 (예시 5개):")
        for f, msg in bad[:5]:
            print("  -", f, "=>", msg)
    return ok

files = filter_valid_netcdf(files)
print(f"[info] 검증 통과 files={len(files)}개")
if not files:
    raise RuntimeError("열 수 있는 파일이 없음.")


# ---------------- [05] 시간축/변수 메타 파악 ----------------
def detect_time_name(nc: NC):
    if ref_time_name and ref_time_name in nc.variables:
        return ref_time_name
    for cand in time_candidates:
        if cand in nc.variables:
            return cand
    raise KeyError(f"시간 좌표 후보({time_candidates}) 중 찾지 못함.")

with NC(files[0], "r") as nc0:
    time_name = detect_time_name(nc0)
    tvar0 = nc0.variables[time_name]
    time_units = getattr(tvar0, "units", "days since 2000-01-01 00:00:00")
    time_cal   = getattr(tvar0, "calendar", "proleptic_gregorian")

    var_meta = {}
    used_dims = set()
    for v in VARS:
        if v not in nc0.variables:
            continue
        var0 = nc0.variables[v]
        dims = var0.dimensions
        if time_name not in dims:
            print(f"[warn] 첫 파일: 변수 {v}에 시간축 {time_name} 없음 → 스킵")
            continue
        t_axis = dims.index(time_name)
        space_dims = tuple(d for d in dims if d != time_name)
        space_shape = tuple(nc0.dimensions[d].size for d in space_dims)
        var_meta[v] = {
            "dims": dims,
            "t_axis": t_axis,
            "space_dims": space_dims,
            "space_shape": space_shape,
        }
        used_dims.update(space_dims)

print(f"[info] monthly mean using time: {time_name}")
print(f"[info] time units={time_units}, calendar={time_cal}")
print(f"[info] 사용 변수 메타: {list(var_meta.keys())}")

if not var_meta:
    raise RuntimeError("월평균 낼 유효 변수가 없음.")


# ---------------- [06] 결과 netCDF 생성 및 스트리밍 누적 ----------------
outfile = outdir / f"surface_monthly_avg_{case}.nc"
print(f"[save init] → {outfile}")

with NC(outfile, "w", format="NETCDF4") as nco:
    # 차원 생성 (time: unlimited)
    nco.createDimension(time_name, None)

    # 공간 차원 및 좌표 변수 복사
    with NC(files[0], "r") as nc0:
        for d in used_dims:
            if d in nc0.dimensions:
                nco.createDimension(d, nc0.dimensions[d].size)

        for d in used_dims:
            if d in nc0.variables:
                v0 = nc0.variables[d]
                vd = nco.createVariable(d, v0.datatype, v0.dimensions)
                vd[:] = v0[:]
                for aname in v0.ncattrs():
                    setattr(vd, aname, getattr(v0, aname))

    # 시간 변수
    tv = nco.createVariable(time_name, "f8", (time_name,))
    tv.units = time_units
    tv.calendar = time_cal
    tv.long_name = "time since initialization"
    tv.field = "time, scalar, series"
    tv.axis = "T"

    # 데이터 변수 생성 (FillValue는 createVariable에서만 설정)
    out_vars = {}
    with NC(files[0], "r") as nc0:
        for v, meta in var_meta.items():
            space_dims = meta["space_dims"]
            out_dims = (time_name,) + space_dims

            v0 = nc0.variables[v]
            fillv = getattr(v0, "_FillValue", None)

            if fillv is not None:
                vo = nco.createVariable(
                    v, "f4", out_dims,
                    zlib=True, complevel=4, shuffle=True,
                    fill_value=fillv,
                )
            else:
                vo = nco.createVariable(
                    v, "f4", out_dims,
                    zlib=True, complevel=4, shuffle=True,
                )

            # attrs 복사 (_FillValue는 이미 위에서 처리했으니 제외)
            for aname in v0.ncattrs():
                if aname == "_FillValue":
                    continue
                setattr(vo, aname, getattr(v0, aname))

            out_vars[v] = vo

    # ---- 월별 누적용 버퍼 (현재 월만) ----
    sum_arr = {}
    cnt_arr = {}
    for v, meta in var_meta.items():
        shape = meta["space_shape"]
        sum_arr[v] = np.zeros(shape, dtype="float64")
        cnt_arr[v] = np.zeros(shape, dtype="int32")

    current_key = [None]   # [ (year, month) ]
    month_index = [0]      # [int]

    def flush_current_month():
        """현재 월 누적값을 평균 내어 파일에 기록하고 버퍼 리셋"""
        if current_key[0] is None:
            return

        y, mth = current_key[0]
        mid = pd.Timestamp(year=y, month=mth, day=15, hour=0, minute=0, second=0)
        tnum = date2num(mid.to_pydatetime(), units=time_units, calendar=time_cal)
        tv[month_index[0]] = tnum

        for v, meta in var_meta.items():
            s = sum_arr[v]
            c = cnt_arr[v]
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = s / c
            mean[c == 0] = np.nan
            out_vars[v][month_index[0], ...] = mean.astype("float32")

            # 다음 달 대비 리셋
            s.fill(0.0)
            c.fill(0)

        print(f"[write] month_index={month_index[0]}, ym={y}-{mth:02d}")
        month_index[0] += 1

    # ---- 파일/시간 순회 (스트리밍) ----
    for path in files:
        print(f"[proc] {path}")
        with NC(path, "r") as nc:
            if time_name not in nc.variables:
                raise KeyError(f"{path}: 시간 변수 {time_name!r} 없음")

            tvar = nc.variables[time_name]
            tvals = tvar[:]
            times = num2date(tvals, units=time_units, calendar=time_cal)

            present_vars = [v for v in VARS if (v in nc.variables and v in var_meta)]
            if not present_vars:
                continue

            ntime = tvar.shape[0]
            for it in range(ntime):
                dt = times[it]
                ym = (dt.year, dt.month)

                # 월 변경 감지
                if current_key[0] is None:
                    current_key[0] = ym
                elif ym != current_key[0]:
                    if ym < current_key[0]:
                        raise RuntimeError(f"시간이 역순으로 감: prev={current_key[0]}, now={ym}")
                    flush_current_month()
                    current_key[0] = ym

                # 현재 step 데이터 누적
                for v in present_vars:
                    meta = var_meta[v]
                    var = nc.variables[v]
                    t_axis = meta["t_axis"]

                    slicer = [slice(None)] * var.ndim
                    slicer[t_axis] = it
                    data = var[tuple(slicer)]  # masked array 가능

                    arr = np.array(data, dtype="float64")
                    fv = getattr(var, "_FillValue", None)
                    if fv is not None:
                        arr[arr == fv] = np.nan
                    mv = getattr(var, "missing_value", None)
                    if mv is not None:
                        arr[arr == mv] = np.nan

                    mask = np.isfinite(arr)
                    s = sum_arr[v]
                    c = cnt_arr[v]
                    s[mask] += arr[mask]
                    c[mask] += 1

    # 마지막 월 flush
    flush_current_month()

print(f"[done] → {outfile}")
