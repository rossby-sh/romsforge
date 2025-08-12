# cmems_ostia_download.py
import os
import yaml
from datetime import datetime as dt
import copernicusmarine

# config 읽기
with open("config_all.yaml", "r") as f:
    cfg = yaml.safe_load(f)

# 경로 준비: base_dir/ostia
base_dir = cfg["output"]["base_dir"].rstrip("/")
target_dir = os.path.join(base_dir, "ostia")
os.makedirs(target_dir, exist_ok=True)
os.chdir(target_dir)  # 여기로 저장되게

# 영역/기간
lat_min, lat_max = cfg["region"]["lat"]
lon_min, lon_max = cfg["region"]["lon"]

start_date = cfg["time"]["start"]  # "YYYY-MM-DD"
end_date   = cfg["time"]["end"]    # "YYYY-MM-DD"
start_iso = f"{start_date}T00:00:00"
end_iso   = f"{end_date}T00:00:00"

# OSTIA L4 NRT SST (METOFFICE) 서브셋
copernicusmarine.subset(
    dataset_id="METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2",
    variables=["analysed_sst", "analysis_error", "mask"],
    minimum_longitude=lon_min,
    maximum_longitude=lon_max,
    minimum_latitude=lat_min,
    maximum_latitude=lat_max,
    start_datetime=start_iso,
    end_datetime=end_iso,
)

print(f"[done] saved into: {target_dir}")

