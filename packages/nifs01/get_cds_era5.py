import os, random, time, yaml
from datetime import datetime as dt, timedelta
import cdsapi
import numpy as np
import multiprocessing as mp

# 워커 (요청 + 리트라이 + 원자적 저장)
def worker(dataset, request, outpath, retries=5):
    client = cdsapi.Client()
    tmp_path = outpath + ".part"
    for i in range(retries):
        try:
            print(f"[{os.getpid()}] retrieve → {outpath} ({request.get('variable')})")
            client.retrieve(dataset, request, tmp_path)
            os.replace(tmp_path, outpath)
            return outpath
        except Exception as e:
            if os.path.exists(tmp_path):
                try: os.remove(tmp_path)
                except: pass
            wait = (2**i) + random.random()
            print(f"[retry {i+1}/{retries}] {e} → sleep {wait:.1f}s")
            time.sleep(wait)
    raise RuntimeError(f"failed: {outpath}")

# 24시간
TIME = [f"{h:02d}:00" for h in np.arange(0,24,3)]

if __name__ == "__main__":
    with open("config_clm.yaml") as f:
        cfg = yaml.safe_load(f)

    start_dt = dt.strptime(cfg["time"]["start"], "%Y-%m-%d")
    end_dt   = dt.strptime(cfg["time"]["end"],   "%Y-%m-%d")
    pad_before = int(cfg["time"].get("pad_before_days", 0))
    pad_after  = int(cfg["time"].get("pad_after_days", 0))

    lat_min, lat_max = cfg["region"]["lat"]
    lon_min, lon_max = cfg["region"]["lon"]
    base_dir = cfg["output"]["base_dir"]

    tmp_dir = os.path.join(base_dir, "era5")
    os.makedirs(tmp_dir, exist_ok=True)

    # 패딩 포함 요청 범위(풀데이로 받을 거라 date는 범위 문자열)
    req_start = (start_dt - timedelta(days=pad_before)).date()
    req_end   = (end_dt   + timedelta(days=pad_after)).date()
    period_tag = f"{req_start.strftime('%Y%m%d')}-{req_end.strftime('%Y%m%d')}"

    # 데이터셋
    dataset_sl = "reanalysis-era5-single-levels"
    dataset_pl = "reanalysis-era5-pressure-levels"
    area = [lat_max, lon_min, lat_min, lon_max]  # N, W, S, E

    # 변수 목록
    accum_vars = [
        "total_precipitation",
        "surface_latent_heat_flux",
        "surface_net_solar_radiation",
        "surface_net_thermal_radiation",
        "surface_solar_radiation_downwards",
        "surface_thermal_radiation_downwards",
        "evaporation",
        "potential_evaporation",
        "runoff",
    ]
    instant_vars = [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_dewpoint_temperature",
        "2m_temperature",
        "mean_sea_level_pressure",
        "sea_surface_temperature",
        "surface_pressure",
        "skin_temperature",
        "total_cloud_cover",
    ]
    q1000_vars = ["specific_humidity"]
    q1000_extra = {"pressure_level": ["1000"]}

    # 공통 요청 빌더 (풀데이, 범위문자열)
    def build_req(var_name, extra=None):
        req = {
            "product_type": ["reanalysis"],
            "variable": [var_name],
            "date": f"{req_start:%Y-%m-%d}/{req_end:%Y-%m-%d}",
            "time": TIME,
            "format": "grib",
            "area": area,
        }
        if extra:
            req.update(extra)
        return req

    jobs = []

    def add_job_var(dataset_name, kind, var_name, extra=None):
        outpath = os.path.join(tmp_dir, f"{kind}__{var_name}_{period_tag}.grib")
        if not os.path.exists(outpath):
            req = build_req(var_name, extra=extra)
            jobs.append((dataset_name, req, outpath))
        else:
            print(f"[skip] exists: {outpath}")

    # single-levels: 누적/순간 변수 각각 파일 하나씩
    for v in accum_vars:
        add_job_var(dataset_sl, "accum", v)
    for v in instant_vars:
        add_job_var(dataset_sl, "inst", v)

    # pressure-levels: q1000
    for v in q1000_vars:
        add_job_var(dataset_pl, "q1000", v, extra=q1000_extra)

    print(f"총 작업 {len(jobs)}개")
    if jobs:
        # 과도한 병렬은 429 잘 남. 기본 3, 환경변수로 조절 가능
        nproc = int(os.getenv("CDS_WORKERS", "1"))
        with mp.Pool(processes=nproc) as pool:
            pool.starmap(worker, jobs, chunksize=1)

    print("--- Done ---")

