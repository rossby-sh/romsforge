import os
import cdsapi
import multiprocessing
import time, random
import yaml
from datetime import datetime as dt, timedelta

def retrieve_with_retry(dataset, request, target, retries=5):
    client = cdsapi.Client()
    for i in range(retries):
        try:
            client.retrieve(dataset, request, target)  # GRIB은 타깃 파일 지정
            return
        except Exception as e:
            wait = (2 ** i) + random.random()
            print(f"[retry {i+1}/{retries}] {e} → sleep {wait:.1f}s")
            time.sleep(wait)
    raise RuntimeError(f"CDS retrieve failed after {retries} retries: {target}")

def build_request_common(years, months, days, area, variables):
    return {
        "product_type": ["reanalysis"],
        "variable": variables,
        "year": years,
        "month": months,
        "day": days,
        "time": ["00:00","03:00","06:00","09:00","12:00","15:00","18:00","21:00"],
        "format": "grib",  # GRIB
        "area": [area[0], area[1], area[2], area[3]]  # [N, W, S, E]
    }

def download_job(args):
    dataset, request, outpath = args
    print(f"--- Download {outpath} ---")
    retrieve_with_retry(dataset, request, outpath)
    return outpath

def generate_time_fields(start_date, end_date):
    current = start_date
    years = set(); months = set(); days = set()
    while current <= end_date:
        years.add(current.strftime("%Y"))
        months.add(current.strftime("%m"))
        days.add(current.strftime("%d"))
        current += timedelta(days=1)
    return sorted(years), sorted(months), sorted(days)

if __name__ == "__main__":
    with open("config.yml") as f:
        cfg = yaml.safe_load(f)

    start_date = dt.strptime(cfg["time"]["start"], "%Y-%m-%d")
    end_date   = dt.strptime(cfg["time"]["end"],   "%Y-%m-%d")
    lat_min, lat_max = cfg["region"]["lat"]
    lon_min, lon_max = cfg["region"]["lon"]
    base_dir = cfg["output"]["base_dir"]

    # 출력 디렉토리 준비
    tmp_dir = os.path.join(base_dir, "tmp_frc")
    os.makedirs(tmp_dir, exist_ok=True)

    years, months, days = generate_time_fields(start_date, end_date)

    print("=== Download ERA5 (GRIB, parallel) ===")
    dataset = "reanalysis-era5-single-levels"

    # 누적계열
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

    # 순간계열
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

    req_accum = build_request_common(
        years, months, days,
        [lat_max, lon_min, lat_min, lon_max],
        accum_vars
    )
    req_instant = build_request_common(
        years, months, days,
        [lat_max, lon_min, lat_min, lon_max],
        instant_vars
    )

    accum_path = os.path.join(tmp_dir, "accum.grib")
    inst_path  = os.path.join(tmp_dir, "inst.grib")

    jobs = [
        (dataset, req_accum, accum_path),
        (dataset, req_instant, inst_path),
    ]

    with multiprocessing.Pool(processes=2) as pool:
        _ = pool.map(download_job, jobs)

    print("--- Done downloads ---")

