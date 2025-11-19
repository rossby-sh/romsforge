
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
surface_raw.py
- 오버랩(4일) 스킵
- 표층 선택
- 월평균 없음
- 시간 원본 그대로 유지
- config의 vars만 추출
- netCDF4 엔진으로 단일 파일 저장
- 처리 과정 로깅 추가
"""

import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import xarray as xr
import dask
from netCDF4 import Dataset as NC

# --- libs/utils.py 로드 ---
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import utils as tl  # noqa: E402


# ---------------- [01] 설정 로드 ----------------
cfg = tl.parse_config("./config_proc_nft.yaml")
_ = tl.load_roms_grid(getattr(cfg, "grdname", None))

case = getattr(cfg, "case", "case")
data_dir = Path(getattr(cfg, "data_dir", ".")).expanduser().resolve()
base_dir = Path(getattr(cfg, "base_dir", ".")).expanduser().resolve()
outdir = base_dir / str(case)
outdir.mkdir(parents=True, exist_ok=True)

data_prefix_header = str(getattr(cfg, "data_prefix_header", ""))

# 📌 ★ 변수는 config에서 가져옴
VARS = list(getattr(cfg, "vars", []))

engine = "netcdf4"
chunks_per_cycle = int(getattr(cfg, "chunks_per_cycle", 8))
overlap = int(getattr(cfg, "overlap", 4))
validate = bool(getattr(cfg, "validate", True))

vn = getattr(cfg, "var_names", None) or {}
time_candidates = list(getattr(vn, "time_candidates", ["ocean_time", "time", "t"]))
lat_candidates  = set(getattr(vn, "lat_candidates", ["lat", "latitude", "y", "eta", "eta_rho"]))
lon_candidates  = set(getattr(vn, "lon_candidates", ["lon", "longitude", "x", "xi", "xi_rho"]))
z_candidates    = list(getattr(vn, "z_candidates", ["s_rho", "z", "depth", "lev", "layer"]))
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

print(f"[info] 전체 파일 {len(files_all)}개 발견")
print(f"[info] 예시 파일: {files_all[:3]}")
print(f"[info] 추출 변수: {VARS}")


# ---------------- [03] 오버랩 제거 ----------------
NAME_RE = re.compile(r"(?P<cycle>\d{4})_(?P<chunk>\d{4})\.nc$")

def select_nonoverlap(all_files, chunks_per_cycle=8, overlap=4):
    groups = {}
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

files = select_nonoverlap(files_all, chunks_per_cycle, overlap)
print(f"[info] 오버랩 제거 후 파일 {len(files)}개")


# ---------------- [04] 유효 파일 필터 ----------------
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
        print(f"[warn] 손상 파일 {len(bad)}개 제외:")
        for f, msg in bad[:5]:
            print("  -", f, "=>", msg)

    return ok

files = filter_valid_netcdf(files)
if not files:
    raise RuntimeError("열 수 있는 파일이 없음.")


# ---------------- [05] 표층 전처리 + 로깅 ----------------
def detect_zdim(da: xr.DataArray):
    for cand in z_candidates:
        if cand in da.dims and da.sizes[cand] > 1:
            return cand
    not_z = set(time_candidates) | lat_candidates | lon_candidates
    for d in da.dims:
        if d not in not_z and da.sizes[d] > 1:
            return d
    return None

def preprocess_surface(ds: xr.Dataset) -> xr.Dataset:
    fname = ds.encoding.get("source", None)
    if fname:
        print(f"[preprocess] {Path(fname).name} 처리 중...")

    keep = [v for v in VARS if v in ds.data_vars]
    if not keep:
        return xr.Dataset()

    out = []
    for v in keep:
        da = ds[v]
        zdim = detect_zdim(da)
        out.append(
            da if zdim is None else da.isel({zdim: -1}, drop=True)
        )
    return xr.merge([x.to_dataset(name=v) for x, v in zip(out, keep)])


# ---------------- [06] 병합 ----------------
print("[info] open_mfdataset 시작...")
ds = xr.open_mfdataset(
    files,
    combine="by_coords",
    preprocess=preprocess_surface,
    parallel=False,
    data_vars="minimal",
    coords="minimal",
    compat="override",
    engine=engine,
)
print("[info] open_mfdataset 완료")


# ---------------- 시간 좌표 찾기 ----------------
time_name = None
if ref_time_name and ref_time_name in ds.coords:
    time_name = ref_time_name
else:
    for cand in time_candidates:
        if cand in ds.coords:
            time_name = cand
            break

if time_name is None:
    raise KeyError("시간 좌표(time 변수)를 찾지 못함.")

ds = ds.sortby(time_name).chunk({time_name: 1})
print(f"[info] 시간 좌표: {time_name}, shape={ds[time_name].shape}")


# ---------------- [07/08] 월평균 없음 → 그대로 사용 ----------------
m = ds


# ---------------- [09] 저장 ----------------
outfile = outdir / f"surface_raw_{case}.nc"
print("[info] NetCDF 쓰기 시작...")

encoding = {v: {"zlib": True, "complevel": 4, "shuffle": True} for v in m.data_vars}

with dask.config.set(scheduler="single-threaded"):
    m.to_netcdf(outfile, encoding=encoding, engine=engine, format="NETCDF4")

print("[done] 저장 완료 →", outfile)
