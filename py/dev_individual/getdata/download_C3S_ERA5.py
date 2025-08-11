import cdsapi
import xarray as xr
import yaml
from datetime import datetime, timedelta

with open("config.yml") as f:
    cfg = yaml.safe_load(f)

print("=== Download ERA5 ===")

# 공통 파라미터
start_date = datetime.strptime(cfg["time"]["start"], "%Y-%m-%d")
end_date   = datetime.strptime(cfg["time"]["end"],   "%Y-%m-%d")
lat_min, lat_max = cfg["region"]["lat"]
lon_min, lon_max = cfg["region"]["lon"]
base_dir = cfg["output"]["base_dir"]

prefix = "ERA5_3hourly"

filename = f"{prefix}_{start_date.strftime('%y%m%d')}-{end_date.strftime('%y%m%d')}.nc"

def generate_time_fields(start_date, end_date):
    current = start_date
    years = set()
    months = set()
    days = set()

    while current <= end_date:
        years.add(current.strftime("%Y"))
        months.add(current.strftime("%m"))
        days.add(current.strftime("%d"))
        current += timedelta(days=1)

    return sorted(years), sorted(months), sorted(days)

years, months, days = generate_time_fields(start_date, end_date)

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "total_precipitation",
        "surface_latent_heat_flux",
        "surface_net_solar_radiation",
        "surface_net_thermal_radiation",
        "surface_sensible_heat_flux",
        "surface_solar_radiation_downwards",
        "surface_thermal_radiation_downwards",
        "evaporation",
        "potential_evaporation",
        "runoff"
    ],
    "year": years,
    "month": months,
    "day": days,
    "time": [
        "00:00", "03:00", "06:00",
        "09:00", "12:00", "15:00",
        "18:00", "21:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [lat_max, lon_min,  # 북,서,남,동
             lat_min, lon_max]
}

print("--- Download accumlates ---")
client = cdsapi.Client()
client.retrieve(dataset, request).download(base_dir+"tmp_frc/accum.nc")

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_dewpoint_temperature",
        "2m_temperature",
        "mean_sea_level_pressure",
        "sea_surface_temperature",
        "surface_pressure",
        "skin_temperature",
        "total_cloud_cover",
    ],
    "year": years,
    "month": months,
    "day": days,
    "time": [
        "00:00", "03:00", "06:00",
        "09:00", "12:00", "15:00",
        "18:00", "21:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [lat_max, lon_min,  # 북,서,남,동
             lat_min, lon_max]
}

print("--- Download instant ---")

client = cdsapi.Client()
client.retrieve(dataset, request).download(base_dir+"tmp_frc/inst.nc")


print("--- Merge and save ERA5 to "+base_dir+filename+" ---")
# Merge files
accm=xr.open_dataset(base_dir+'tmp_frc/accum.nc').loc[dict(latitude=slice(60,5),longitude=slice(100,170))] 
inst=xr.open_dataset(base_dir+'tmp_frc/inst.nc').loc[dict(latitude=slice(60,5),longitude=slice(100,170))] 

ERA5 = xr.merge([accm,inst],compat='override')
ERA5.to_netcdf(base_dir+filename)














