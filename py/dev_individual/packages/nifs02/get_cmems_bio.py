# cmems_bio_download_mp.py
from datetime import datetime
import os, yaml, multiprocessing as mp
import copernicusmarine  # pip install copernicusmarine

def download_one(args):
    item, lon, lat, start_iso, end_iso, depth_min, depth_max, out_dir, period_tag = args
    out_name = f"CMEMS_{item['name']}_{period_tag}.nc"
    copernicusmarine.subset(
        dataset_id=item["dataset_id"],
        variables=item["variables"],
        minimum_longitude=lon[0],
        maximum_longitude=lon[1],
        minimum_latitude=lat[0],
        maximum_latitude=lat[1],
        start_datetime=start_iso,
        end_datetime=end_iso,
        minimum_depth=depth_min,
        maximum_depth=depth_max,
        output_filename=out_name,
        output_directory=out_dir,
        force_download=True,
    )
    return f"[OK] {item['name']}: {os.path.join(out_dir, out_name)}"

if __name__ == "__main__":
    # config 읽기
    with open("config_all.yaml") as f:
        cfg = yaml.safe_load(f)

    # 공통 설정
    start_date = datetime.strptime(cfg["time"]["start"], "%Y-%m-%d")
    end_date   = datetime.strptime(cfg["time"]["end"],   "%Y-%m-%d")
    latitude   = tuple(cfg["region"]["lat"])   # (min_lat, max_lat)
    longitude  = tuple(cfg["region"]["lon"])   # (min_lon, max_lon)

    depth_min = 0.0
    depth_max = 5902.05810546875

    # 출력 디렉토리: base_dir/cmems_bio
    base_dir   = cfg["output"]["base_dir"].rstrip("/")
    output_dir = os.path.join(base_dir, "cmems_bio")
    os.makedirs(output_dir, exist_ok=True)

    # 기간(ISO)
    start_iso = f"{start_date:%Y-%m-%d}T00:00:00"
    end_iso   = f"{end_date:%Y-%m-%d}T00:00:00"
    period_tag = f"{start_date:%Y%m%d}-{end_date:%Y%m%d}"

    # 항목 정의 (각 항목당 파일 1개)
    download_configs = [
        {"name": "car", "dataset_id": "cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m", "variables": ["dissic", "ph", "talk"]},
        {"name": "co2", "dataset_id": "cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m", "variables": ["spco2"]},
        {"name": "nut", "dataset_id": "cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m", "variables": ["fe", "no3", "po4", "si"]},
        {"name": "pft", "dataset_id": "cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m", "variables": ["chl", "phyc"]},
        {"name": "bio", "dataset_id": "cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m", "variables": ["o2", "nppv"]},
    ]

    # 풀 작업 인자 구성
    jobs = [
        (item, longitude, latitude, start_iso, end_iso, depth_min, depth_max, output_dir, period_tag)
        for item in download_configs
    ]

    # 병렬 실행 (2~3 권장)
    procs = int(os.environ.get("CMEMS_WORKERS", "3"))
    with mp.Pool(processes=procs) as pool:
        for msg in pool.imap_unordered(download_one, jobs, chunksize=1):
            print(msg)

