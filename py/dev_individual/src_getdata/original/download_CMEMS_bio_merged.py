from datetime import datetime
from dateutil.relativedelta import relativedelta
import copernicusmarine

# 공통 설정
start_date = datetime(2025, 7, 1)
end_date = datetime(2025, 7, 30)
delta_t = relativedelta(months=1)

longitude = (100, 170)
latitude = (5, 60)
depth_min = 0.5057600140571594
depth_max = 5902.05810546875

# 다운로드 항목 설정 리스트
download_configs = [
    {
        "name": "car",
        "wpth": "/data/share/DATA/RAW/Bvar/Carbon/",
        "dataset_id": "cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m",
        "variables": ["dissic", "ph", "talk"]
    },
    {
        "name": "co2",
        "wpth": "/data/share/DATA/RAW/Bvar/CO2/",
        "dataset_id": "cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m",
        "variables": ["spco2"]
    },
    {
        "name": "nut",
        "wpth": "/data/share/DATA/RAW/Bvar/NUT/",
        "dataset_id": "cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m",
        "variables": ["fe", "no3", "po4", "si"]
    },
    {
        "name": "pft",
        "wpth": "/data/share/DATA/RAW/Bvar/PFT/",
        "dataset_id": "cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m",
        "variables": ["chl", "phyc"]
    },
    {
        "name": "bio",
        "wpth": "/data/share/DATA/RAW/Bvar/BIO/",
        "dataset_id": "cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m",
        "variables": ["o2", "nppv"]
    }
]

# 다운로드 루프
def download_monthly_data(config):
    current_date = start_date
    while current_date <= end_date:
        start_datetime = current_date
        end_datetime = current_date + delta_t

        output_filename = f"CMEMS_data_{config['name']}_{start_datetime.strftime('%Y-%m')}.nc"

        copernicusmarine.subset(
            dataset_id=config["dataset_id"],
            variables=config["variables"],
            minimum_longitude=longitude[0],
            maximum_longitude=longitude[1],
            minimum_latitude=latitude[0],
            maximum_latitude=latitude[1],
            start_datetime=start_datetime,
            end_datetime=end_datetime,
            minimum_depth=depth_min,
            maximum_depth=depth_max,
            output_filename=output_filename,
            output_directory=config["wpth"],
            force_download=True
        )

        print(f"[X] {config['name'].upper()}: {start_datetime.date()} ~ {end_datetime.date()} → {output_filename}")

        current_date += delta_t

# 실행
for config in download_configs:
    download_monthly_data(config)

