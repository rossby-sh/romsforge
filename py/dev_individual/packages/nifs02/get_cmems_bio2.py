
# get_cmems_bio.py (안전 패치)
from datetime import datetime
import os, yaml, time
import copernicusmarine

def download_one(item, lon, lat, start_iso, end_iso, out_dir, period_tag, max_retry=6):
    out_name = f"CMEMS_{item['name']}_{period_tag}.nc"
    wait = 2.0
    for i in range(max_retry):
        try:
            copernicusmarine.subset(
                dataset_id=item["dataset_id"],
                variables=item["variables"],
                minimum_longitude=lon[0], maximum_longitude=lon[1],
                minimum_latitude=lat[0],  maximum_latitude=lat[1],
                start_datetime=start_iso, end_datetime=end_iso,
                # depth는 제품마다 축 다름 → 일단 제거 (필요하면 항목별로 따로 지정)
                output_filename=out_name,
                output_directory=out_dir,
            )
            return f"[OK] {item['name']}: {os.path.join(out_dir, out_name)}"
        except Exception as e:
            msg = str(e)
            # 429나 S3 rate 관련이면 백오프 재시도
            if "Too Many Requests" in msg or "429" in msg:
                time.sleep(wait)
                wait = min(wait * 2, 60.0)
                continue
            return f"[ERR] {item['name']}: {msg}"
    return f"[ERR] {item['name']}: repeated 429 / rate-limited"

if __name__ == "__main__":
    with open("config_all.yaml") as f:
        cfg = yaml.safe_load(f)

    start_date = datetime.strptime(cfg["time"]["start"], "%Y-%m-%d")
    end_date   = datetime.strptime(cfg["time"]["end"],   "%Y-%m-%d")

    # 끝시각을 하루 끝으로 (빈 구간 방지)
    start_iso = f"{start_date:%Y-%m-%d}T00:00:00"
    end_iso   = f"{end_date:%Y-%m-%d}T23:59:59"
    period_tag = f"{start_date:%Y%m%d}-{end_date:%Y%m%d}"

    latitude   = tuple(cfg["region"]["lat"])
    longitude  = tuple(cfg["region"]["lon"])

    base_dir   = cfg["output"]["base_dir"].rstrip("/")
    output_dir = os.path.join(base_dir, "cmems_bio")
    os.makedirs(output_dir, exist_ok=True)

    # (선택) 환경변수 로그인
    USER = os.getenv("COPERNICUSMARINE_USERNAME")
    PWD  = os.getenv("COPERNICUSMARINE_PASSWORD")
    if USER and PWD:
        copernicusmarine.login(USER, PWD)

    download_configs = [
        {"name": "car", "dataset_id": "cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m", "variables": ["dissic", "ph", "talk"]},
        {"name": "co2", "dataset_id": "cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m", "variables": ["spco2"]},
        {"name": "nut", "dataset_id": "cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m", "variables": ["fe", "no3", "po4", "si"]},
        {"name": "pft", "dataset_id": "cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m", "variables": ["chl", "phyc"]},
        {"name": "bio", "dataset_id": "cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m", "variables": ["o2", "nppv"]},
    ]

    # 완전 직렬 실행 + 항목 사이에 약간 쉼 주기(과요청 방지)
    for item in download_configs:
        print(download_one(item, longitude, latitude, start_iso, end_iso, output_dir, period_tag))
        time.sleep(1.0)
