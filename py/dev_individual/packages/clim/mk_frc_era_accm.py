
# =================================
# grib_accq_region.py  ACCUM+Q 전용 (subset 지원)
# =================================
# - 입력: accum__*.grib, q1000__*.grib
# - 출력 변수: srf, lwrad_down, rain, Qair (+ time, lat, lon)
# - 저장은 하지 않음. 아래 forcing_vars_accq 만 제공

import os, sys, glob, yaml, numpy as np, eccodes as ec
from datetime import datetime
from datetime import datetime as dt

# libs
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
from utils import compute_relative_time

with open('config_all.yaml', 'r') as f:
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
    AREA = (
        float(a.get("north")),
        float(a.get("west")),
        float(a.get("south")),
        float(a.get("east")),
    )

# 재사용 유틸 (간단 복사)

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
    lat_mask = (lat >= min(S,N)) & (lat <= max(S,N))
    lat_idx = np.where(lat_mask)[0]
    if lat_idx.size == 0:
        raise RuntimeError('area: 위도 범위 격자 없음')
    lat_slice = slice(lat_idx[0], lat_idx[-1]+1)
    if W <= E:
        lon_mask = (lon >= W) & (lon <= E)
        lon_idx = np.where(lon_mask)[0]
        if lon_idx.size == 0:
            raise RuntimeError('area: 경도 범위 격자 없음')
        lon_slices = [slice(lon_idx[0], lon_idx[-1]+1)]
        lon_out = lon[lon_slices[0]]
    else:
        idx1 = np.where(lon >= W)[0]
        idx2 = np.where(lon <= E)[0]
        lon_slices = []
        if idx1.size:
            lon_slices.append(slice(idx1[0], idx1[-1]+1))
        if idx2.size:
            lon_slices.append(slice(idx2[0], idx2[-1]+1))
        lon_out = np.concatenate([lon[s] for s in lon_slices])
    return lat_slice, lon_slices, lat[lat_slice], lon_out


def msg_time(gid):
    vdate = ec.codes_get(gid, 'validityDate')
    vtime = ec.codes_get(gid, 'validityTime')
    hh = int(vtime) // 100
    return np.datetime64(datetime(int(str(vdate)[:4]), int(str(vdate)[4:6]), int(str(vdate)[6:8]), hh), 's')


def step_hours(gid):
    hours = 1
    if ec.codes_is_defined(gid, 'stepRange'):
        sr = str(ec.codes_get(gid, 'stepRange'))
        try:
            a, b = sr.split('-'); hours = max(1, int(b) - int(a))
        except Exception:
            try:
                hours = max(1, int(sr))
            except Exception:
                hours = 1
    return hours

# 파일 스캔
p_acc = sorted(glob.glob(os.path.join(ERA5_DIR, 'accum__*.grib')))
p_q   = sorted(glob.glob(os.path.join(ERA5_DIR, 'q1000__*.grib')))
if not p_acc:
    raise SystemExit('accum__*.grib 없음')
if not p_q:
    print('[warn] q1000__*.grib 없음 (Qair 미출력)')

# anchor: 강수 파일 우선, 없으면 다른 accum
anchor = None
for p in p_acc:
    if 'total_precipitation' in os.path.basename(p):
        anchor = p; break
if anchor is None:
    anchor = p_acc[0]

LAT_FULL, LON_FULL, LAT_REV = read_grid_from_first_msg(anchor)
LAT_SL, LON_SLICES, LAT, LON = subset_indices(LAT_FULL, LON_FULL, AREA)

# 시간축 (anchor)
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

epoch = anchor_times.astype('datetime64[s]').astype('int64')
TIME = compute_relative_time(epoch, 'seconds since 1970-01-01 00:00:00', MY_TIME_REF)

# 누적 변환: 에너지 J m-2 → W m-2, 물기둥 m → m s-1
ENERGY = {'ssr','ssrd','str','strd','slhf','sshf'}
MASS   = {'tp','e','pev','ro'}

srf = lwrad_down = rain = None
Qair = None

# accum 처리
for path in p_acc:
    with open(path, 'rb') as f:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None: continue
        short = ec.codes_get(gid, 'shortName'); ec.codes_release(gid)

    frames = []
    with open(path, 'rb') as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None: break
            try:
                t = msg_time(gid)
                hours = step_hours(gid)
                arr = np.asarray(ec.codes_get_values(gid), dtype=np.float64)
                Nj = ec.codes_get(gid, 'Nj'); Ni = ec.codes_get(gid, 'Ni')
                arr = arr.reshape(Nj, Ni)
                lat1= ec.codes_get(gid, 'latitudeOfFirstGridPointInDegrees')
                latN= ec.codes_get(gid, 'latitudeOfLastGridPointInDegrees')
                if lat1 > latN: arr = arr[::-1, :]
                if len(LON_SLICES) == 1:
                    arr = arr[LAT_SL, LON_SLICES[0]]
                else:
                    arr = np.concatenate([arr[LAT_SL, sl] for sl in LON_SLICES], axis=1)
                # 변환
                if short in ENERGY:
                    arr = arr / (hours * 3600.0)
                elif short in MASS:
                    arr = arr / (hours * 3600.0)
                frames.append((t, arr))
            finally:
                ec.codes_release(gid)
    if not frames:
        continue
    frames.sort(key=lambda x: x[0])
    data = np.stack([a for _, a in frames], axis=0)

    if short == 'ssr':
        srf = data
    elif short == 'strd':
        lwrad_down = data
    elif short == 'tp':
        rain = data

# q1000 처리(Qair kg/kg)
for path in p_q:
    with open(path, 'rb') as f:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None: continue
        short = ec.codes_get(gid, 'shortName'); ec.codes_release(gid)
    if short != 'q':
        continue
    frames = []
    with open(path, 'rb') as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None: break
            try:
                t = msg_time(gid)
                arr = np.asarray(ec.codes_get_values(gid), dtype=np.float64)
                Nj = ec.codes_get(gid, 'Nj'); Ni = ec.codes_get(gid, 'Ni')
                arr = arr.reshape(Nj, Ni)
                lat1= ec.codes_get(gid, 'latitudeOfFirstGridPointInDegrees')
                latN= ec.codes_get(gid, 'latitudeOfLastGridPointInDegrees')
                if lat1 > latN: arr = arr[::-1, :]
                if len(LON_SLICES) == 1:
                    arr = arr[LAT_SL, LON_SLICES[0]]
                else:
                    arr = np.concatenate([arr[LAT_SL, sl] for sl in LON_SLICES], axis=1)
                frames.append((t, arr))
            finally:
                ec.codes_release(gid)
    if not frames:
        continue
    frames.sort(key=lambda x: x[0])
    data = np.stack([a for _, a in frames], axis=0)
    Qair = data

# 결과 변수명 매핑
srf_values = srf
lwrad_down_values = lwrad_down
rain_values = rain  # m s-1
qair_values = Qair  # kg kg-1

forcing_vars_accq = {
    'Qair':        qair_values,
    'srf':         srf_values,
    'lwrad_down':  lwrad_down_values,
    'rain':        rain_values,
}

# 여기서 저장하지 않음. TIME/LON/LAT/forcing_vars_accq 를 외부 저장 코드로 넘겨 사용.
