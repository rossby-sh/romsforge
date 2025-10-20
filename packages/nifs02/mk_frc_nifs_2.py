
# grib_min.py  GRIB 처리(저수준 eccodes) + valid_time 정렬 + NetCDF 저장
import os
import sys
import yaml
import numpy as np
import eccodes as ec
from datetime import datetime, timedelta
from datetime import datetime as dt
from netCDF4 import num2date

# libs
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'libs')))
import create_F as cn
from utils import compute_relative_time


# ======================
# 설정 읽기
# ======================
with open('config_all.yaml', 'r') as f:
    config = yaml.safe_load(f)

base = config["output"]["base_dir"].rstrip("/")
start = dt.strptime(config["time"]["start"], "%Y-%m-%d")
end   = dt.strptime(config["time"]["end"],   "%Y-%m-%d")
period_tag = f"{start:%Y%m%d}-{end:%Y%m%d}"

os.makedirs(base + "/roms_inputs", exist_ok=True)

INST_FILE  = os.path.join(base, "era5", f"inst_{period_tag}.grib")
ACCUM_FILE = os.path.join(base, "era5", f"accum_{period_tag}.grib")
Q1000_FILE = os.path.join(base, "era5", f"q1000_{period_tag}.grib")  # 없으면 넘어감
W_FILE     = os.path.join(base, "roms_inputs", f"nifs_frc_{period_tag}.nc")
MY_TIME_REF = config["time_ref"]  # 예: "hours since 2000-01-01 00:00:00"

def must_exist(path, name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"[{name}] 없음: {path}")

must_exist(INST_FILE,  "inst")
must_exist(ACCUM_FILE, "accum")
if not os.path.exists(Q1000_FILE):
    Q1000_FILE = None

print("[inst ]", INST_FILE)
print("[accum]", ACCUM_FILE)
print("[q1000]", Q1000_FILE or "(none)")
print("[out  ]", W_FILE)


# ======================
# 유틸: (t, arr) 리스트 → 배열 패킹
# ======================
def _pack_series(series):
    """series: {shortName: [(t, arr), ...]} → {shortName: {'times': np.datetime64[s], 'data': (T,Ny,Nx)}}"""
    out = {}
    for k, lst in series.items():
        lst.sort(key=lambda x: x[0])  # 시간 정렬
        ts  = np.array([t for t, _ in lst], dtype='datetime64[s]')
        dat = np.stack([a for _, a in lst], axis=0)  # (T, Ny, Nx)

        # 중복 time 제거(첫 등장 유지)
        uniq, idx = np.unique(ts, return_index=True)
        if len(uniq) != len(ts):
            keep = np.sort(idx)
            ts = ts[keep]
            dat = dat[keep, ...]
        out[k] = {"times": ts, "data": dat}
    return out


# ======================
# GRIB 리더 (한 파일 전체)
# - instant: validityDate/Time == 관측 시각(time)
# - accum : validityDate/Time == 유효 시각(valid_time = base time + step)
#           → 여기서 만든 t는 '끝 시각'이므로 나중 슬라이스에서 그대로 사용
# - accum=True일 때는 stepRange로 시간구간 길이(hours) 추출해 레이트로 환산
# ======================
ACC_M_TO_MM   = {"e", "pev", "ro"}                        # m over N hours → mm/hr
ACC_J_TO_WM2  = {"ssr", "ssrd", "str", "strd", "slhf", "sshf"}  # J/m^2 over N hours → W/m^2

def read_grib(path, wanted=None, accum=False):
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
                    lat_rev = (lat1 > latN)  # 보통 ERA5는 북→남 스캔 → True

                # 유효 시각(끝 시각; accum은 valid_time, instant는 관측시각)
                vdate = ec.codes_get(gid, "validityDate")   # 예) 20250901
                vtime = ec.codes_get(gid, "validityTime")   # 예) 0, 300, 2100
                hh = int(vtime) // 100
                t  = np.datetime64(datetime(int(str(vdate)[:4]),
                                            int(str(vdate)[4:6]),
                                            int(str(vdate)[6:8]), hh), "s")

                # 필드
                arr = np.asarray(ec.codes_get_values(gid), dtype=np.float64).reshape(Nj, Ni)

                # 위도 방향 정렬 통일: 남→북(증가) 되도록 뒤집기
                if lat_rev:
                    arr = arr[::-1, :]

                # 누적 → 레이트 환산 (메시지별 stepRange 사용)
                if accum and (short in ACC_M_TO_MM or short in ACC_J_TO_WM2):
                    hours = 1
                    if ec.codes_is_defined(gid, "stepRange"):
                        sr = str(ec.codes_get(gid, "stepRange"))
                        try:
                            # "0-1", "1-3" 같은 형식
                            a, b = sr.split("-")
                            hours = max(1, int(b) - int(a))
                        except Exception:
                            # "1" 같은 단일 값
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

    # 좌표 벡터 생성(요청 영역은 0.25° 격자 가정)
    lon = lon1 + dlon * np.arange(Ni)
    if lat_rev:
        # arr를 뒤집었으니 위도도 남→북(작은→큰) 순으로 맞춰줌
        lat_start = min(lat1, latN)
        lat = lat_start + abs(dlat) * np.arange(Nj)
    else:
        lat = lat1 + dlat * np.arange(Nj)

    vars_dict = _pack_series(series)
    return lat, lon, vars_dict


# ======================
# 데이터 읽기
# ======================
# instant
inst_wanted = {"10u","10v","2t","2d","msl","sst","skt","sp","tcc"}
inst_lat, inst_lon, INS = read_grib(INST_FILE, wanted=inst_wanted, accum=False)

# accum (누적 → 레이트 환산 포함; 유효시간은 끝 시각으로 t에 저장되어 있음)
accum_wanted = {"e","pev","ro","slhf","ssr","ssrd","str","strd","tp"}
acc_lat, acc_lon, ACC = read_grib(ACCUM_FILE, wanted=accum_wanted, accum=True)

# q1000 (선택)
Q = {}
if Q1000_FILE:
    if os.path.exists(Q1000_FILE) and os.path.getsize(Q1000_FILE) > 0:
        _, _, Q = read_grib(Q1000_FILE, wanted={"q"}, accum=False)

# 좌표(instant 기준 사용)
LON = inst_lon
LAT = inst_lat


# ======================
# 시간창 슬라이스 & 공통 시간축 정렬
# - instant:  start 00:00 ~ end 23:00
# - accum  :  start 01:00 ~ (end+1d) 00:00  (누적의 끝 시각이므로 +1h 쉬프트됨)
# ======================
def _slice_by_time(ts, arr, t0, t1):
    m = (ts >= t0) & (ts <= t1)
    return ts[m], arr[m]

def _align_to(ref_ts, ts, arr):
    pos = {t: i for i, t in enumerate(ts)}
    sel_idx = [pos[t] for t in ref_ts if t in pos]
    keep_ts = np.array([t for t in ref_ts if t in pos], dtype='datetime64[s]')
    if len(sel_idx) == 0:
        # 비어 있으면 빈 배열 반환 (호출부에서 체크)
        return keep_ts, arr[:0, ...]
    return keep_ts, arr[sel_idx, ...]

# 타깃 기간 정의
t0_inst = np.datetime64(f"{start:%Y-%m-%d}T00:00")
t1_inst = np.datetime64(f"{end:%Y-%m-%d}T23:00")

t0_acc  = np.datetime64(f"{start:%Y-%m-%d}T01:00")
t1_acc  = np.datetime64(f"{(end + timedelta(days=1)):%Y-%m-%d}T00:00")

# 기준축: instant의 2t
ts_inst_ref = INS["2t"]["times"].copy()
ts_inst_ref, _ = _slice_by_time(ts_inst_ref, INS["2t"]["data"], t0_inst, t1_inst)

# 모든 instant 변수: 슬라이스 → 기준축으로 리인덱스
for k in list(INS.keys()):
    ts, da = INS[k]["times"], INS[k]["data"]
    ts, da = _slice_by_time(ts, da, t0_inst, t1_inst)
    ts_new, da_new = _align_to(ts_inst_ref, ts, da)
    INS[k]["times"], INS[k]["data"] = ts_new, da_new

# accum의 가능한 시간축(예: tp 기준 하나 꺼냄) → 누적윈도우로 슬라이스
if len(ACC) == 0:
    raise RuntimeError("ACC(누적) 변수를 찾지 못했어. accum GRIB을 확인해줘.")
sample_k = "tp" if "tp" in ACC else next(iter(ACC.keys()))
ts_acc_any = ACC[sample_k]["times"]
ts_acc_any, _ = _slice_by_time(ts_acc_any, ACC[sample_k]["data"], t0_acc, t1_acc)

# 공통 시간축(= instant 기준축 ∩ accum-윈도우)
common_ts = np.array(sorted(set(ts_inst_ref) & set(ts_acc_any)), dtype='datetime64[s]')
if common_ts.size == 0:
    raise RuntimeError("공통 시간 교집합이 비었어. 요청 기간/버퍼/자르기 규칙을 확인해줘.")

# instant를 공통축으로 재정렬
for k in list(INS.keys()):
    ts, da = INS[k]["times"], INS[k]["data"]
    _, da = _align_to(common_ts, ts, da)
    INS[k]["times"], INS[k]["data"] = common_ts, da

# accum을 공통축으로 재정렬 (accum은 이미 t=유효시각이므로 그냥 교집합 맞추기)
for k in list(ACC.keys()):
    ts, da = ACC[k]["times"], ACC[k]["data"]
    ts, da = _slice_by_time(ts, da, t0_acc, t1_acc)
    _, da = _align_to(common_ts, ts, da)
    ACC[k]["times"], ACC[k]["data"] = common_ts, da

# q1000도 instant 창에 맞춰 정렬(있을 때)
if "q" in Q:
    ts, da = Q["q"]["times"], Q["q"]["data"]
    ts, da = _slice_by_time(ts, da, t0_inst, t1_inst)
    _, da = _align_to(common_ts, ts, da)
    Q["q"]["times"], Q["q"]["data"] = common_ts, da

print("=== 공통 시간축 ===")
print(common_ts[0], "→", common_ts[-1], "N=", common_ts.size)


# ======================
# 변수 매핑 (네가 쓰던 이름 유지)
# ======================
def _get(name, src):
    if name not in src:
        raise KeyError(f"'{name}' 변수를 찾지 못했어.")
    return src[name]["data"]

# instant
u_values   = _get("10u", INS)
v_values   = _get("10v", INS)
t2k_values = _get("2t",  INS)                 # K
t2_values  = t2k_values - 273.15              # °C
d2_values  = _get("2d",  INS) - 273.15        # °C
pair_values= _get("msl", INS) / 100.0         # Pa → hPa
cloud_values = _get("tcc", INS)
sst_values = _get("sst", INS) - 273.15        # °C (후처리)
sst_values[sst_values > 1000] = 0.0           # 이상치 클립
TK_values  = t2k_values

# accum (이미 레이트 환산 끝난 상태)
srf_values        = _get("ssr",  ACC)         # W m-2
lwrad_down_values = _get("strd", ACC)         # W m-2
lwrad_values      = _get("str",  ACC)         # W m-2
# tp는 read_grib에서 mm/hr로 환산 안 하고 W/m2만 환산했으므로 여기서 강수 레이트로 만들자:
# 실제 ACC에서 tp는 'm over N hours'임 → 위에서 hours로 나눠주지 않으면 평균강수강도 안 됨.
# 위 read_grib에서 ACC_M_TO_MM에 "tp"가 없었음 → 추가 변환
if "tp" in ACC:
    # ACC["tp"]["data"]는 'm over hours'에서 (환산 안 됐다면) 그대로일 수 있으므로 확인
    # 안전하게 stepRange 분해 없이 시간차분 대신, 이미 t=valid_time으로 맞췄으니
    # 여기선 m/s로 바꿔서 ROMS용 kg m-2 s-1로 쓰자: mm/hr → m/s는 (mm/hr)/1000/3600
    # 만약 read_grib에서 ACC_M_TO_MM에 'tp' 추가해서 mm/hr로 바꿨다면 아래 줄만 사용하면 됨.
    rain_values = ACC["tp"]["data"] * 1000.0 / 3600.0  # m over 1h → mm/hr → m/s
else:
    raise KeyError("'tp' (total_precipitation)가 ACC에서 누락됐어.")

# qair @ 1000hPa (선택)
if Q1000_FILE and "q" in Q:
    qair_values = Q["q"]["data"] * 1000.0  # kg/kg → g/kg (혹은 모델 요구 포맷에 맞춰 조정)
else:
    raise RuntimeError("q1000(grib) 파일에서 'q'를 찾지 못했어. config_era5.yml의 input_q1000를 확인해줘.")

# RH (상대습도, %)
rh_values = 100.0 * (
    np.exp((17.625 * d2_values) / (243.04 + d2_values)) /
    np.exp((17.625 * t2_values) / (243.04 + t2_values))
)

# dqdsst (필요하면 주석 해제)
wind_speed = np.sqrt(u_values**2 + v_values**2)
# dqdsst_values = cn.get_dqdsst(sst_values, TK_values, 1.2, wind_speed, qair_values / 1000.0)

# 모두 동일한 시간축인지 체크
T_expected = common_ts.size
for name, arr in {
    "Uwind": u_values, "Vwind": v_values, "Tair": t2_values, "Qair": qair_values,
    "RH": rh_values, "Cloud": cloud_values, "sst": sst_values, "srf": srf_values,
    "lwrad": lwrad_values, "lwrad_down": lwrad_down_values, "rain": rain_values,
    "Pair": pair_values,
}.items():
    assert arr.shape[0] == T_expected, f"{name} 시간길이 불일치: {arr.shape[0]} vs {T_expected}"

# ======================
# 시간 변환: 공통 시간축 → ref_time
# ======================
times64 = common_ts.astype("datetime64[s]")
epoch_sec = times64.astype("int64")  # seconds since 1970-01-01
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
    # 'dqdsst':      dqdsst_values,
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
