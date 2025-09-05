# modis_aqua_chl_nrt_download_mp.py
from datetime import datetime, timedelta
import os, yaml, multiprocessing as mp
import xarray as xr

def dates_between(start: datetime, end: datetime):
    d = start
    while d <= end:
        yield d
        d += timedelta(days=1)

def worker(args):
    one_date, out_dir, lon_min, lon_max, lat_min, lat_max = args
    yyyy = one_date.strftime('%Y')
    mmdd = one_date.strftime('%m%d')
    ymd  = one_date.strftime('%Y%m%d')

    url = (
        f"https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/"
        f"{yyyy}/{mmdd}/AQUA_MODIS.{ymd}.L3m.DAY.CHL.chlor_a.9km.NRT.nc"
    )
    out_file = os.path.join(out_dir, f"AQUA_MODIS.{ymd}.NW_Pacific.chlor_a_NRT.nc")

    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return f"[SKIP] {ymd} (exists)"

    try:
        # OPeNDAP 원격 열기 → 영역 잘라 저장
        with xr.open_dataset(url) as ds:
            chl = ds["chlor_a"].sel(
                lon=slice(lon_min, lon_max),
                lat=slice(lat_max, lat_min)  # 위도는 보통 남북 반대방향으로 슬라이스
            )
            tmp = out_file + ".part"
            chl.to_netcdf(tmp)
            os.replace(tmp, out_file)
        return f"[OK]  {ymd} → {out_file}"
    except Exception as e:
        return f"[FAIL] {ymd}: {e}"

if __name__ == "__main__":
    # config 읽기
    with open("config_all.yaml") as f:
        cfg = yaml.safe_load(f)

    start_date = datetime.strptime(cfg["time"]["start"], "%Y-%m-%d")
    end_date   = datetime.strptime(cfg["time"]["end"],   "%Y-%m-%d")
    lat_min, lat_max = cfg["region"]["lat"]   # [5, 60]
    lon_min, lon_max = cfg["region"]["lon"]   # [100, 180]

    base_dir = cfg["output"]["base_dir"].rstrip("/")
    out_dir  = os.path.join(base_dir, "modis")
    os.makedirs(out_dir, exist_ok=True)

    procs = int(os.environ.get("MODIS_WORKERS", "3"))

    jobs = [(d, out_dir, lon_min, lon_max, lat_min, lat_max)
            for d in dates_between(start_date, end_date)]

    print(f"총 {len(jobs)}일, 프로세스 {procs}개 → {out_dir}")

    with mp.Pool(processes=procs) as pool:
        for msg in pool.imap_unordered(worker, jobs, chunksize=1):
            print(msg)

