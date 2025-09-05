
# grib_min.py  one-layer parser (eccodes + numpy)
import eccodes as ec
import numpy as np
from datetime import datetime

import sys
import os
import datetime as dt
import yaml
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_F as cn
from utils import compute_relative_time

# Load config
with open('config_era5.yml', 'r') as f:
    config = yaml.safe_load(f)

# Construct file paths
R_FILE = os.path.join(config['input_file'])
W_FILE = os.path.join(config['output_file'])



# 누적 에너지(J m-2) → 평균 플럭스(W m-2)로 바꿀 변수들
ENERGY_ACCUM = {"ssr", "strd", "slhf", "sshf"}

def read_grib(path, wanted=None, accum=False):
    """
    GRIB 파일을 읽어 변수별 시계열을 반환.
    - wanted: {'u10','v10','t2m',...} 지정하면 그 변수만 읽음(없으면 전부)
    - accum=True면 누적 변수를 dt로 나눠 구간량/평균 플럭스로 변환:
        * ENERGY_ACCUM → W m-2
        * tp → mm/hr
        * 그 외 누적 → 원 단위의 구간량(diff)
    return: lat(Ny), lon(Nx), vars:{shortName: {'times':datetime64[T],'data':(T,Ny,Nx)}}
    """
    series = {}
    Ni = Nj = None
    lat1 = latN = lon1 = dlat = dlon = None
    lat_rev = False

    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                short = ec.codes_get(gid, "shortName")
                if wanted and short not in wanted:
                    continue
                if Ni is None:
                    Ni  = ec.codes_get(gid, "Ni")
                    Nj  = ec.codes_get(gid, "Nj")
                    lat1= ec.codes_get(gid, "latitudeOfFirstGridPointInDegrees")
                    latN= ec.codes_get(gid, "latitudeOfLastGridPointInDegrees")
                    lon1= ec.codes_get(gid, "longitudeOfFirstGridPointInDegrees")
                    dlon= ec.codes_get(gid, "iDirectionIncrementInDegrees")
                    dlat= ec.codes_get(gid, "jDirectionIncrementInDegrees")
                    lat_rev = (lat1 > latN)  # ERA5 보통 북→남

                vdate = ec.codes_get(gid, "validityDate")  # 20250701
                vtime = ec.codes_get(gid, "validityTime")  # 0, 300, 2100
                hh = int(vtime) // 100
                t = datetime(int(str(vdate)[:4]), int(str(vdate)[4:6]), int(str(vdate)[6:8]), hh)

                vals = ec.codes_get_values(gid)
                arr  = np.asarray(vals, dtype=np.float64).reshape(Nj, Ni)
                if lat_rev:  # 위도 오름차순(남→북) 정렬
                    arr = arr[::-1, :]
                series.setdefault(short, []).append((t, arr))
            finally:
                ec.codes_release(gid)

    if Ni is None:
        raise RuntimeError(f"No messages in {path}")

    # 좌표
    lon = (lon1 + np.arange(Ni) * dlon) % 360.0
    if lat_rev:
        lat_top, lat_bot = latN, lat1
    else:
        lat_top, lat_bot = lat1, latN
    lat = np.linspace(lat_top, lat_bot, Nj)

    # 변수별 스택 및 (필요 시) 누적 → 구간 변환
    out = {}
    for short, lst in series.items():
        lst.sort(key=lambda x: x[0])
        times = np.array([np.datetime64(x[0], "s") for x in lst])
        data  = np.stack([x[1] for x in lst], axis=0)  # (T, Ny, Nx)

        if accum and data.shape[0] >= 2:
            tsec = times.astype("int64")  # epoch seconds
            dt   = np.diff(tsec)          # (T-1,)
            inc  = np.diff(data, axis=0)  # (T-1,Ny,Nx)

            if short == "tp":
                # m 누적 → mm/hr
                inc = inc * 1000.0 / (dt[:, None, None] / 3600.0)
            elif short in ENERGY_ACCUM:
                # J m-2 누적 → W m-2
                inc = inc / dt[:, None, None]
            # 첫 타임 패드(NaN)
            pad = np.full((1,)+inc.shape[1:], np.nan, inc.dtype)
            data = np.concatenate([pad, inc], axis=0)

        out[short] = {"times": times, "data": data}

    return lat, lon, out


pth='/home/shjo/ROMS/romsforge/py/dev_individual/getdata/test/tmp_frc/'

accum_lat, accum_lon, ACC = read_grib(pth+"accum_20250630-20250801.grib", accum=True)
inst_lat,  inst_lon,  INS = read_grib(pth+"inst_20250630-20250801.grib",  accum=False)
_, _, Q1000 = read_grib(pth+"q1000_20250630-20250801.grib", wanted={"q"})  # specific_humidity shortName이 'q'임

'''
def list_grib_vars(path):
    names=set()
    with open(path,'rb') as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None: break
            try:
                names.add(ec.codes_get(gid, "shortName"))
            finally:
                ec.codes_release(gid)
    return sorted(names)

print(list_grib_vars(pth+"inst_20250630-20250801.grib"))
print(list_grib_vars(pth+"accum_20250630-20250801.grib"))
print(list_grib_vars(pth+"q1000_20250630-20250801.grib"))
'''

# 예시: 바람/기온/압력
sst_values = INS["sst"]["data"]          # (T,Ny,Nx)
u_values = INS["10u"]["data"]          # (T,Ny,Nx)
v_values = INS["10v"]["data"]
t2_values = INS["2t"]["data"] - 273.15
d2_values = INS["2d"]["data"] - 273.15
pair_values = INS["msl"]["data"] / 100.0  # Pa → hPa
cloud_values = INS["tcc"]["data"]
# 예시: 복사/강수 (accum → 이미 변환됨)
srf_values  = ACC["ssr"]["data"]     # W m-2
lwrad_down_values = ACC["strd"]["data"]    # W m-2
lwrad_values = ACC["str"]["data"]    # W m-2
rain_values = ACC["tp"]["data"]      # mm hr-1

# q(1000hPa)
qair_values = Q1000["q"]["data"]        # kg/kg at 1000hPa

TK_values        =  INS["2t"]["data"]

# Calculate Qair
rh_values = 100 * (
    np.exp((17.625 * d2_values) / (243.04 + d2_values)) /
    np.exp((17.625 * t2_values) / (243.04 + t2_values))
)


# Calculate dqdsst
wind_speed = np.sqrt(u_values ** 2 + v_values ** 2)
dqdsst_values = cn.get_dqdsst(
    sst_values,
    TK_values,
    1.2,
    wind_speed,
    qair_values / 1000
)
MY_TIME_REF = config['ref_time']
# Convert time
TIME_CONVERTED_NUM = compute_relative_time(TIMES[:], TIMES.units, MY_TIME_REF)

print("=== Convert time ===")
print(num2date(TIME_CONVERTED_NUM[:5],MY_TIME_REF))
print(num2date(TIME_CONVERTED_NUM[-5:],MY_TIME_REF))

# Assemble output variables
forcing_vars = {
    'Uwind':       u_values,
    'Vwind':       v_values,
    'Tair':        t2_values,
    'Qair':        qair_values,
    'RH'   :       rh_values,
    'Cloud':       cloud_values,
    'sst':         sst_values,
    'dqdsst':      dqdsst_values,
    'srf':         srf_values,
    'lwrad':       lrf_values,
    'lwrad_down':  lrf_down_values,
    'rain':        rain_values,
    'Pair':        pair_values,
}

print("=== Save to netcdf ===")

# Write output NetCDF
cn.createF_era5(
    W_FILE,
    LON,
    LAT,
    TIME_CONVERTED_NUM,
    MY_TIME_REF,
    forcing_vars,
    "NETCDF3_64BIT_OFFSET"
)









