# grib_min_multi.py — 변수별 GRIB 파일 입력(여러 개) 지원 버전
# - 기존: inst.accum.q1000 3개 큰 GRIB 파일 → 변경: 변수별 GRIB 파일 각각 존재
# - 동작: era5 폴더 내 파일명을 스캔하여 유형(inst/accum/q1000)별로 자동 수집
# - 조건: q1000의 'q'가 없으면 종료(추정 금지)

import os
import sys
import glob
import yaml
import numpy as np
import eccodes as ec
from datetime import datetime
from netCDF4 import num2date
from datetime import datetime as dt

# libs
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_F as cn
from utils import compute_relative_time

# ======================
# 설정 읽기
# ======================
with open('config_single.yaml', 'r') as f:
    config = yaml.safe_load(f)

base = config["output"]["base_dir"].rstrip("/")
start = dt.strptime(config["time"]["start"], "%Y-%m-%d")
end   = dt.strptime(config["time"]["end"],   "%Y-%m-%d")
period_tag = f"{start:%Y%m%d}-{end:%Y%m%d}"

os.makedirs(base+"/roms_inputs", exist_ok=True)
ERA5_DIR = os.path.join(base, "era5")
W_FILE   = os.path.join(base, "roms_inputs", f"nifs_frc_{period_tag}.nc")
MY_TIME_REF = config["time_ref"]  # 예: "hours since 2000-01-01 00:00:00"

print("[scan ]", ERA5_DIR)

# ======================
# 파일 스캐너 (변수별 파일)
# ======================
# 기대 파일명 패턴 예시:
#   inst__2m_temperature_19900101-20250731.grib
#   accum__total_precipitation_19900101-20250731.grib
#   q1000__specific_humidity_19900101-20250731.grib
# 접두어(inst__/accum__/q1000__)로 유형을 판별하고, 각 파일의 shortName을 키로 맵 구성

def scan_variable_files(folder):
    paths = glob.glob(os.path.join(folder, "*.grib"))
    inst_files = {}
    accum_files = {}
    q1000_files = {}
    for p in sorted(paths):
        name = os.path.basename(p)
        if name.startswith("inst__"):
            inst_files[name] = p
        elif name.startswith("accum__"):
            accum_files[name] = p
        elif name.startswith("q1000__"):
            q1000_files[name] = p
    return inst_files, accum_files, q1000_files

# ======================
# GRIB 리더 (단일 파일)
# ======================
# - 파일 내 단일 shortName만 있다고 가정(변수별 파일이므로)
# - accum=True면 메시지(stepRange) 길이로 나눠서 율로 환산(기존 로직 유지)

#ACC_M_TO_MM = {"e", "pev", "ro"}
#ACC_J_TO_WM2 = {"ssr", "ssrd", "str", "strd", "slhf", "sshf"}
ACC_M_TO_MM = {"tp"}
ACC_J_TO_WM2 = {"ssr", "strd"}



def read_grib_single(path, accum=False):
    series = {}
    Ni = Nj = None
    lat1 = latN = lon1 = dlat = dlon = None
    lat_rev = False
    short_name = None

    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                short = ec.codes_get(gid, "shortName")
                if short_name is None:
                    short_name = short

                if Ni is None:
                    Ni  = ec.codes_get(gid, "Ni")
                    Nj  = ec.codes_get(gid, "Nj")
                    lat1= ec.codes_get(gid, "latitudeOfFirstGridPointInDegrees")
                    latN= ec.codes_get(gid, "latitudeOfLastGridPointInDegrees")
                    lon1= ec.codes_get(gid, "longitudeOfFirstGridPointInDegrees")
                    dlon= ec.codes_get(gid, "iDirectionIncrementInDegrees")
                    dlat= ec.codes_get(gid, "jDirectionIncrementInDegrees")
                    lat_rev = (lat1 > latN)

                vdate = ec.codes_get(gid, "validityDate")
                vtime = ec.codes_get(gid, "validityTime")
                hh = int(vtime) // 100
                t  = np.datetime64(datetime(int(str(vdate)[:4]),
                                            int(str(vdate)[4:6]),
                                            int(str(vdate)[6:8]), hh), "s")

                arr = np.asarray(ec.codes_get_values(gid), dtype=np.float64).reshape(Nj, Ni)
                if lat_rev:
                    arr = arr[::-1, :]

                if accum and (short in ACC_M_TO_MM or short in ACC_J_TO_WM2):
                    hours = 1
                    if ec.codes_is_defined(gid, "stepRange"):
                        sr = str(ec.codes_get(gid, "stepRange"))
                        try:
                            a, b = sr.split("-")
                            hours = max(1, int(b) - int(a))
                        except Exception:
                            try:
                                hours = max(1, int(sr))
                            except Exception:
                                hours = 1
                    if short in ACC_M_TO_MM:
                        arr = arr * 1000.0 / hours          # m over N hours → mm/hr
                    elif short in ACC_J_TO_WM2:
                        arr = arr / (hours * 3600.0)        # J m-2 over N hours → W m-2

                series.setdefault(short, []).append((t, arr))
            finally:
                if 'gid' in locals() and gid is not None:
                    ec.codes_release(gid)

    if Ni is None:
        raise RuntimeError(f"No messages in {path}")

    lon = (lon1 + np.arange(Ni) * dlon) % 360.0
    lat_top, lat_bot = (latN, lat1) if lat_rev else (lat1, latN)
    lat = np.linspace(lat_top, lat_bot, Nj)

    # 단일 shortName만 있다고 전제
    if len(series) != 1:
        raise RuntimeError(f"{path}: shortName가 하나가 아닙니다: {list(series.keys())}")

    short = list(series.keys())[0]
    lst = series[short]
    lst.sort(key=lambda x: x[0])
    times = np.array([t for t,_ in lst], dtype="datetime64[s]")
    data  = np.stack([a for _,a in lst], axis=0)

    return lat, lon, short, {short: {"times": times, "data": data}}

# ======================
# 스캔 & 로드
# ======================
inst_files, accum_files, q1000_files = scan_variable_files(ERA5_DIR)

if not inst_files:
    raise RuntimeError("inst__* GRIB 파일을 찾지 못했습니다.")
if not accum_files:
    raise RuntimeError("accum__* GRIB 파일을 찾지 못했습니다.")
if not q1000_files:
    print("[warn ] q1000__* 파일이 없습니다. 'q'가 필요하다면 이후 단계에서 실패합니다.")

INS = {}
ACC = {}
Q   = {}
LON = LAT = None

# instant
for name, path in inst_files.items():
    lat, lon, short, d = read_grib_single(path, accum=False)
    if LON is None:
        LON, LAT = lon, lat
    else:
        if not(np.allclose(LON, lon) and np.allclose(LAT, lat)):
            raise RuntimeError(f"격자 불일치(inst): {name}")
    INS.update(d)
    print(f"[inst ] {short:>5} ← {name}")

# accum
for name, path in accum_files.items():
    lat, lon, short, d = read_grib_single(path, accum=True)
    if not(np.allclose(LON, lon) and np.allclose(LAT, lat)):
        raise RuntimeError(f"격자 불일치(accum): {name}")
    ACC.update(d)
    print(f"[accum] {short:>5} ← {name}")

# q1000 (선택)
for name, path in q1000_files.items():
    lat, lon, short, d = read_grib_single(path, accum=False)
    if short != 'q':
        print(f"[warn ] q1000 파일이지만 shortName={short} (무시하지 않고 로드) ← {name}")
    if not(np.allclose(LON, lon) and np.allclose(LAT, lat)):
        raise RuntimeError(f"격자 불일치(q1000): {name}")
    Q.update(d)
    print(f"[q1000] {short:>5} ← {name}")

# 필수 변수 존재 확인(기존 스키마 유지)
required_ins = {"10u","10v","2t","msl","tcc"}
missing_ins = sorted(list(required_ins - set(INS.keys())))
if missing_ins:
    raise RuntimeError(f"instant 누락: {missing_ins}")

required_acc = {"ssr","strd","tp"}  # 최소 세트(필요시 확장)
missing_acc = sorted(list(required_acc - set(ACC.keys())))
if missing_acc:
    raise RuntimeError(f"accum 누락: {missing_acc}")

if 'q' not in Q:
    raise RuntimeError("q1000(grib) 파일에서 'q'를 찾지 못했습니다. (추정 금지)")

# ======================
# 변수 매핑 (동일 로직 유지)
# ======================
#sst_values = INS["sst"]["data"]                         # K → 아래서 변환
u_values   = INS["10u"]["data"]
v_values   = INS["10v"]["data"]
t2_values  = INS["2t"]["data"]  - 273.15               # °C
#d2_values  = INS["2d"]["data"]  - 273.15               # °C
#TK_values  = INS["2t"]["data"]                          # K
pair_values= INS["msl"]["data"] / 100.0                # Pa → hPa
cloud_values = INS["tcc"]["data"]

srf_values       = ACC["ssr"]["data"]                  # W m-2 (accum→율 환산 완료)
lwrad_down_values= ACC["strd"]["data"]                 # W m-2
#lwrad_values     = ACC["str"]["data"]                  # W m-2
rain_values      = ACC["tp"]["data"] * 1000 / 3600     # mm hr-1 → m/s(=kg m-2 s-1)

qair_values = Q["q"]["data"] * 1000.0                  # kg/kg → g/kg (기존 코드 호환)

# NaN 처리 등 필요시
#sst_values = sst_values - 273.15                         # °C
#sst_values[sst_values > 1000] = 0.0

# RH
#rh_values = 100.0 * (
#    np.exp((17.625 * d2_values) / (243.04 + d2_values)) /
#    np.exp((17.625 * t2_values) / (243.04 + t2_values))
#)

# dqdsst
#wind_speed = np.sqrt(u_values**2 + v_values**2)
#dqdsst_values = cn.get_dqdsst(
#    sst_values,
#    TK_values,
#    1.2,
#    wind_speed,
#    qair_values / 1000.0   # g/kg → kg/kg
#)

# ======================
# 시간 변환: GRIB → ref_time
# ======================
# 기준 시간: instant의 2t 시각 사용

# 모든 시계열이 같은 시간축을 공유한다고 가정(ERA5 6-hourly)
# 일치 검사는 최소한 2t/10u 정도만 확인

times64 = INS["2t"]["times"]
if not np.array_equal(times64, INS["10u"]["times"]):
    raise RuntimeError("time axis mismatch: 2t vs 10u")

epoch_sec = times64.astype("datetime64[s]").astype("int64")
epoch_units = "seconds since 1970-01-01 00:00:00"
TIME_CONVERTED_NUM = compute_relative_time(epoch_sec, epoch_units, MY_TIME_REF)

print("=== Convert time ===")
print(num2date(TIME_CONVERTED_NUM[:5],  MY_TIME_REF))
print(num2date(TIME_CONVERTED_NUM[-5:], MY_TIME_REF))

# ======================
# 저장
# ======================
forcing_vars = {
    'Uwind':       u_values,
    'Vwind':       v_values,
    'Tair':        t2_values,
    'Qair':        qair_values,
#    'RH':          rh_values,
    'Cloud':       cloud_values,
#    'sst':         sst_values,
#    'dqdsst':      dqdsst_values,
    'srf':         srf_values,
#    'lwrad':       lwrad_values,
    'lwrad_down':  lwrad_down_values,
    'rain':        rain_values,
    'Pair':        pair_values,
}

print("=== Save to netcdf ===")
cn.createF_era5(
    W_FILE,
    LON,
    LAT,
    TIME_CONVERTED_NUM,
    MY_TIME_REF,
    forcing_vars,
    "NETCDF3_64BIT_OFFSET"
)
print("[done]", W_FILE)

