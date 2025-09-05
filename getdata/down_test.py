import os, random, time, yaml
from datetime import datetime as dt, timedelta
import cdsapi
import multiprocessing as mp

# 워커 함수
def worker(dataset, request, outpath, retries=5):
    client = cdsapi.Client()
    for i in range(retries):
        try:
            print(f"[{os.getpid()}] retrieve → {outpath}")
            client.retrieve(dataset, request, outpath)
            return outpath
        except Exception as e:
            wait = (2**i) + random.random()
            print(f"[retry {i+1}/{retries}] {e} → sleep {wait:.1f}s")
            time.sleep(wait)
    raise RuntimeError(f"failed: {outpath}")

# 날짜별 가변 시간 생성
time8 = ["00:00","03:00","06:00","09:00","12:00","15:00","18:00","21:00"]
def hours_for_date(d, start_dt, end_dt):
    hours = time8[:]
    if d.date() == start_dt.date():
        hours = [h for h in hours if int(h[:2]) >= start_dt.hour]
    if d.date() == end_dt.date():
        hours = [h for h in hours if int(h[:2]) <= end_dt.hour]
    return hours

if __name__ == "__main__":
    # config 읽기
    with open("config.yml") as f:
        cfg = yaml.safe_load(f)

    # 시(hour)까지 읽기 (예: "2025-06-30 12")
    start_dt = dt.strptime(cfg["time"]["start"], "%Y-%m-%d %H")
    end_dt   = dt.strptime(cfg["time"]["end"],   "%Y-%m-%d %H")
    lat_min, lat_max = cfg["region"]["lat"]
    lon_min, lon_max = cfg["region"]["lon"]
    base_dir = cfg["output"]["base_dir"]

    tmp_dir = os.path.join(base_dir, "tmp_frc")
    os.makedirs(tmp_dir, exist_ok=True)

    dataset = "reanalysis-era5-single-levels"
    area = [lat_max, lon_min, lat_min, lon_max]

    # 변수 리스트
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

    # 날짜별로 잡 생성
    jobs = []
    cur = start_dt
    while cur.date() <= end_dt.date():
        day_hours = hours_for_date(cur, start_dt, end_dt)
        if not day_hours:
            cur += timedelta(days=1)
            continue

        y = cur.strftime("%Y")
        m = cur.strftime("%m")
        d = cur.strftime("%d")

        # 누적계열
        for var in accum_vars:
            req = {
                "product_type": ["reanalysis"],
                "variable": [var],
                "year": [y], "month": [m], "day": [d],
                "time": day_hours,
                "format": "grib",
                "area": area,
            }
            out = os.path.join(tmp_dir, f"accum_{var}_{y}{m}{d}.grib")
            if not os.path.exists(out):
                jobs.append((dataset, req, out))

        # 순간계열
        for var in instant_vars:
            req = {
                "product_type": ["reanalysis"],
                "variable": [var],
                "year": [y], "month": [m], "day": [d],
                "time": day_hours,
                "format": "grib",
                "area": area,
            }
            out = os.path.join(tmp_dir, f"inst_{var}_{y}{m}{d}.grib")
            if not os.path.exists(out):
                jobs.append((dataset, req, out))

        cur += timedelta(days=1)

    print(f"총 작업 {len(jobs)}개")
    with mp.Pool(processes=2) as pool:
        pool.starmap(worker, jobs, chunksize=1)
    print("--- Done ---")

