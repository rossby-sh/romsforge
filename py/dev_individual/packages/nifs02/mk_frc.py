# grib_min.py — GRIB 처리(저수준 eccodes) + ref_time 변환 + NetCDF 저장
import os
import sys
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
with open('config_frc.yaml', 'r') as f:
    config = yaml.safe_load(f)

base = config["output"]["base_dir"].rstrip("/")
start = dt.strptime(config["time"]["start"], "%Y-%m-%d")
end   = dt.strptime(config["time"]["end"],   "%Y-%m-%d")
period_tag = f"{start:%Y%m%d}-{end:%Y%m%d}"

os.makedirs(base+"/roms_inputs",exist_ok=True)

INST_FILE  = os.path.join(base, "era5", f"inst_{period_tag}.grib")
ACCUM_FILE = os.path.join(base, "era5", f"accum_{period_tag}.grib")
Q1000_FILE = os.path.join(base, "era5", f"q1000_{period_tag}.grib")  # 없으면 넘어가자
W_FILE     = os.path.join(base, "roms_inputs", f"nifs_frc_{period_tag}.nc")
MY_TIME_REF = config["time_ref"]  # 예: "hours since 2000-01-01 00:00:00"

def must_exist(path, name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"[{name}] 없음: {path}")

must_exist(INST_FILE,  "inst")
must_exist(ACCUM_FILE, "accum")
# q1000은 선택
if not os.path.exists(Q1000_FILE):
    Q1000_FILE = None

print("[inst ]", INST_FILE)
print("[accum]", ACCUM_FILE)
print("[q1000]", Q1000_FILE or "(none)")
print("[out  ]", W_FILE)


# ======================
# GRIB 리더 (한 레이어)
# ======================
ENERGY_ACCUM = {"ssr", "ssrd", "str", "strd", "slhf", "sshf"}  # J m-2 누적 → W m-2

def read_grib(path, wanted=None, accum=False):
    """
    GRIB 파일 → 변수별 시계열(dict).
    wanted: {'10u','10v','2t','2d','msl','sst','tcc', ...}
    accum=True: 누적 필드(tp/ssr/str/strd/slhf/sshf/e/pev/ro)를
                메시지별 누적구간 길이(stepRange)로 나눠서
                mm/hr 또는 W/m^2로 환산 (diff 사용 안 함)
    return: lat(Ny), lon(Nx), vars: {name: {'times': datetime64[T], 'data': (T,Ny,Nx)}}
    """
    series = {}
    Ni = Nj = None
    lat1 = latN = lon1 = dlat = dlon = None
    lat_rev = False

    ACC_M_TO_MM = { "e", "pev", "ro"}
    ACC_J_TO_WM2 = {"ssr", "ssrd", "str", "strd", "slhf", "sshf"}

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
                    lat_rev = (lat1 > latN)  # 보통 북→남 스캔

                # 유효 시각
                vdate = ec.codes_get(gid, "validityDate")   # 20250701
                vtime = ec.codes_get(gid, "validityTime")   # 0, 300, 2100
                hh = int(vtime) // 100
                t  = np.datetime64(datetime(int(str(vdate)[:4]),
                                            int(str(vdate)[4:6]),
                                            int(str(vdate)[6:8]), hh), "s")

                # 필드
                arr = np.asarray(ec.codes_get_values(gid), dtype=np.float64).reshape(Nj, Ni)
                if lat_rev:
                    arr = arr[::-1, :]  # 남→북 정렬

                # 누적 환산(메시지 단위로 처리)
                if accum and (short in ACC_M_TO_MM or short in ACC_J_TO_WM2):
                    # 누적 구간 길이(시간) 추출: stepRange="0-1","1-3" 등
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
                        # m over N hours → mm/hr
                        arr = arr * 1000.0 / hours
                    elif short in ACC_J_TO_WM2:
                        # J m-2 over N hours → W m-2
                        arr = arr / (hours * 3600.0)

                series.setdefault(short, []).append((t, arr))

            finally:
                if 'gid' in locals() and gid is not None:
                    ec.codes_release(gid)

    if Ni is None:
        raise RuntimeError(f"No messages in {path}")

    # 좌표
    lon = (lon1 + np.arange(Ni) * dlon) % 360.0
    lat_top, lat_bot = (latN, lat1) if lat_rev else (lat1, latN)
    lat = np.linspace(lat_top, lat_bot, Nj)

    # 시간 정렬 & 스택
    out = {}
    for short, lst in series.items():
        lst.sort(key=lambda x: x[0])
        times = np.array([t for t,_ in lst], dtype="datetime64[s]")
        data  = np.stack([a for _,a in lst], axis=0)  # (T,Ny,Nx)
        out[short] = {"times": times, "data": data}

    return lat, lon, out


# ======================
# 데이터 읽기
# ======================
# instant
inst_wanted = {"10u","10v","2t","2d","msl","skt","sp","sst","tcc"}
inst_lat, inst_lon, INS = read_grib(INST_FILE, wanted=inst_wanted, accum=False)

# accum (누적 → 구간 변환 포함)
accum_wanted = {"e","pev","ro","slhf","ssr","ssrd","str","strd","tp"}
acc_lat, acc_lon, ACC = read_grib(ACCUM_FILE, wanted=accum_wanted, accum=True)

# q1000 (선택)
Q = {}
if Q1000_FILE:
    if os.path.exists(Q1000_FILE) and os.path.getsize(Q1000_FILE) > 0:
        _, _, Q = read_grib(Q1000_FILE, wanted={"q"}, accum=False)

# 좌표 (instant 기준으로 사용)
LON = inst_lon
LAT  = inst_lat

# ======================
# 변수 매핑 (네가 쓰던 이름 유지)
# ======================
sst_values = INS["sst"]["data"]                         # K → 아래서 변환
u_values   = INS["10u"]["data"]
v_values   = INS["10v"]["data"]
t2_values  = INS["2t"]["data"]  - 273.15               # °C
d2_values  = INS["2d"]["data"]  - 273.15               # °C
TK_values  = INS["2t"]["data"]                          # K
pair_values= INS["msl"]["data"] / 100.0                # Pa → hPa
cloud_values = INS["tcc"]["data"]

srf_values       = ACC["ssr"]["data"]                  # W m-2
lwrad_down_values= ACC["strd"]["data"]                 # W m-2
lwrad_values     = ACC["str"]["data"]                  # W m-2
rain_values      = ACC["tp"]["data"] * 1000 / 3600     # mm hr-1 → m/s(=kg m-2 s-1)

qair_values = Q["q"]["data"]*1000 if "q" in Q else None     # g/kg @ 1000hPa (없으면 None)

# NaN 처리 등 필요시
sst_values = sst_values - 273.15                       # °C
sst_values[sst_values>1000] = 0.0

# RH
rh_values = 100.0 * (
    np.exp((17.625 * d2_values) / (243.04 + d2_values)) /
    np.exp((17.625 * t2_values) / (243.04 + t2_values))
)

# dqdsst
wind_speed = np.sqrt(u_values**2 + v_values**2)
if qair_values is None:
    # qair 없으면 근사 필요하지만, 추정 금지 조건이라 그대로 둔다.
    raise RuntimeError("q1000(grib) 파일에서 'q'를 찾지 못했어. config_era5.yml의 input_q1000를 확인해줘.")
dqdsst_values = cn.get_dqdsst(
    sst_values,
    TK_values,
    1.2,
    wind_speed,
    qair_values / 1000.0
)

# ======================
# 시간 변환: GRIB → ref_time
# ======================
# 기준 시간 배열 하나만 사용(instant의 2t 기준)
times64 = INS["2t"]["times"]                     # datetime64[s]
epoch_sec = times64.astype("datetime64[s]").astype("int64")  # seconds since 1970-01-01
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
    'RH':          rh_values,
    'Cloud':       cloud_values,
    'sst':         sst_values,
    'dqdsst':      dqdsst_values,
    'srf':         srf_values,
    'lwrad':       lwrad_values,
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

