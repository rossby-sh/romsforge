#!/usr/bin/env python3
# grib_inst_only.py — ERA5 instant 변수만 처리 (위경도 subset + createF_era5 호출)
# 입력: ERA5 변수별 GRIB (inst__*.grib)
# 출력: ROMS forcing NetCDF (Uwind, Vwind, Tair, Cloud, Pair)
# 필요: eccodes, utils.compute_relative_time, create_F.createF_era5

import os
import sys
import glob
import yaml
import numpy as np
import eccodes as ec
from datetime import datetime
from datetime import datetime as dt

# libs
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
from utils import compute_relative_time
import create_F as cn  # createF_era5

# -----------------
# 설정
# -----------------
with open('config_single.yaml', 'r') as f:
    config = yaml.safe_load(f)

base = config["output"]["base_dir"].rstrip("/")
ERA5_DIR = os.path.join(base, "era5")
start = dt.strptime(config["time"]["start"], "%Y-%m-%d")
end   = dt.strptime(config["time"]["end"],   "%Y-%m-%d")
period_tag = f"{start:%Y%m%d}-{end:%Y%m%d}"
MY_TIME_REF = config["time_ref"]

AREA = None
if "area" in config and config["area"]:
    a = config["area"]
    AREA = (float(a.get("north")), float(a.get("west")),
            float(a.get("south")), float(a.get("east")))

W_FILE = os.path.join(base, "roms_inputs", f"nifs_frc_inst_{period_tag}.nc")
os.makedirs(os.path.dirname(W_FILE), exist_ok=True)

# -----------------
# 유틸
# -----------------
def read_grid_from_first_msg(path):
    with open(path, 'rb') as f:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None:
            raise RuntimeError(f"Empty GRIB: {path}")
        Ni  = ec.codes_get(gid, 'Ni');  Nj  = ec.codes_get(gid, 'Nj')
        lat1= ec.codes_get(gid, 'latitudeOfFirstGridPointInDegrees')
        latN= ec.codes_get(gid, 'latitudeOfLastGridPointInDegrees')
        lon1= ec.codes_get(gid, 'longitudeOfFirstGridPointInDegrees')
        dlon= ec.codes_get(gid, 'iDirectionIncrementInDegrees')
        lat_rev = (lat1 > latN)
        ec.codes_release(gid)
    lon = (lon1 + np.arange(Ni) * dlon) % 360.0
    lat = np.linspace(latN, lat1, Nj) if lat_rev else np.linspace(lat1, latN, Nj)
    return lat, lon, lat_rev

def subset_indices(lat, lon, area):
    if not area:
        return slice(None), [slice(None)], lat, lon
    N, W, S, E = area
    lat_mask = (lat >= min(S, N)) & (lat <= max(S, N))
    lat_idx = np.where(lat_mask)[0]
    if lat_idx.size == 0:
        raise RuntimeError('area: 위도 범위 격자 없음')
    lat_slice = slice(lat_idx[0], lat_idx[-1] + 1)
    if W <= E:
        lon_mask = (lon >= W) & (lon <= E)
        lon_idx = np.where(lon_mask)[0]
        if lon_idx.size == 0:
            raise RuntimeError('area: 경도 범위 격자 없음')
        lon_slices = [slice(lon_idx[0], lon_idx[-1] + 1)]
        lon_out = lon[lon_slices[0]]
    else:
        idx1 = np.where(lon >= W)[0]
        idx2 = np.where(lon <= E)[0]
        lon_slices = []
        if idx1.size:
            lon_slices.append(slice(idx1[0], idx1[-1] + 1))
        if idx2.size:
            lon_slices.append(slice(idx2[0], idx2[-1] + 1))
        lon_out = np.concatenate([lon[s] for s in lon_slices])
    return lat_slice, lon_slices, lat[lat_slice], lon_out

def msg_time(gid):
    vdate = ec.codes_get(gid, 'validityDate')
    vtime = ec.codes_get(gid, 'validityTime')
    hh = int(vtime) // 100
    return np.datetime64(datetime(int(str(vdate)[:4]),
                                  int(str(vdate)[4:6]),
                                  int(str(vdate)[6:8]), hh), 's')

# -----------------
# 파일 스캔
# -----------------
paths = sorted(glob.glob(os.path.join(ERA5_DIR, 'inst__*.grib')))
if not paths:
    raise SystemExit('inst__*.grib 없음')

# anchor: 2m_temperature 우선, 없으면 10m_u
anchor = None
for p in paths:
    if '2m_temperature' in os.path.basename(p):
        anchor = p; break
if anchor is None:
    for p in paths:
        if '10m_u_component_of_wind' in os.path.basename(p):
            anchor = p; break
if anchor is None:
    raise SystemExit('anchor(inst) 없음: 2m_temperature 또는 10m_u 필요')

# 좌표 및 subset 인덱스
LAT_FULL, LON_FULL, LAT_REV = read_grid_from_first_msg(anchor)
LAT_SL, LON_SLICES, LAT, LON = subset_indices(LAT_FULL, LON_FULL, AREA)
Ny, Nx = LAT.size, LON.size

# 시간축 (anchor에서 추출)
anchor_times = []
with open(anchor, 'rb') as f:
    while True:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None: break
        try:
            anchor_times.append(msg_time(gid))
        finally:
            ec.codes_release(gid)
anchor_times = np.array(anchor_times, dtype='datetime64[s]')
anchor_times.sort()
T = anchor_times.size

# epoch → TIME_CONVERTED_NUM
epoch = anchor_times.astype('datetime64[s]').astype('int64')
TIME_CONVERTED_NUM = compute_relative_time(
    epoch,
    'seconds since 1970-01-01 00:00:00',
    MY_TIME_REF
)

# 시간 매핑 dict
time_index = { int(t.astype('datetime64[s]').astype('int64')): i for i, t in enumerate(anchor_times) }

# -----------------
# 필요 변수 목록과 사전할당
# -----------------
need_short = {
    '10u': 'Uwind',
    '10v': 'Vwind',
    '2t' : 'Tair',
    'msl': 'Pair',
    'tcc': 'Cloud',
}

u_values      = np.full((T, Ny, Nx), np.nan, dtype=np.float32)
v_values      = np.full((T, Ny, Nx), np.nan, dtype=np.float32)
t2_values     = np.full((T, Ny, Nx), np.nan, dtype=np.float32)  # degC
pair_values   = np.full((T, Ny, Nx), np.nan, dtype=np.float32)  # hPa
cloud_values  = np.full((T, Ny, Nx), np.nan, dtype=np.float32)  # percent (0..100)

found = set()

# -----------------
# 변수별 파일 처리 (사전할당에 직접 채움)
# -----------------
for path in paths:
    # shortName 확인
    with open(path, 'rb') as f:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None: continue
        short = ec.codes_get(gid, 'shortName')
        ec.codes_release(gid)

    if short not in need_short:
        continue
    found.add(short)

    with open(path, 'rb') as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None: break
            try:
                t = msg_time(gid)
                tkey = int(t.astype('datetime64[s]').astype('int64'))
                if tkey not in time_index:
                    raise RuntimeError(f"time mismatch in {os.path.basename(path)} at {t}")

                idx = time_index[tkey]
                arr = np.asarray(ec.codes_get_values(gid), dtype=np.float32)
                Nj = ec.codes_get(gid, 'Nj'); Ni = ec.codes_get(gid, 'Ni')
                arr = arr.reshape(Nj, Ni)
                # 위도 뒤집기
                lat1= ec.codes_get(gid, 'latitudeOfFirstGridPointInDegrees')
                latN= ec.codes_get(gid, 'latitudeOfLastGridPointInDegrees')
                if lat1 > latN:
                    arr = arr[::-1, :]
                # subset
                if len(LON_SLICES) == 1:
                    arr = arr[LAT_SL, LON_SLICES[0]]
                else:
                    arr = np.concatenate([arr[LAT_SL, sl] for sl in LON_SLICES], axis=1)

                # 변환 및 기록
                if short == '10u':
                    u_values[idx, :, :] = arr  # m/s
                elif short == '10v':
                    v_values[idx, :, :] = arr  # m/s
                elif short == '2t':
                    t2_values[idx, :, :] = arr - 273.15  # K → °C
                elif short == 'msl':
                    pair_values[idx, :, :] = arr / 100.0  # Pa → hPa
                elif short == 'tcc':
                    cloud_values[idx, :, :] = arr * 100.0  # fraction → percent
            finally:
                ec.codes_release(gid)

missing = sorted(list(set(need_short.keys()) - found))
if missing:
    raise RuntimeError(f"instant 누락: {missing}")

# -----------------
# 저장 호출
# -----------------
forcing_vars = {
    'Uwind': u_values,
    'Vwind': v_values,
    'Tair':  t2_values,
    'Cloud': cloud_values,
    'Pair':  pair_values,
}

print("=== Save to netcdf ===")
cn.createF_era5_n3(
    W_FILE,
    LON,
    LAT,
    TIME_CONVERTED_NUM,
    MY_TIME_REF,
    forcing_vars,
    "NETCDF3_64BIT_OFFSET"
)
print("[done]", W_FILE)

