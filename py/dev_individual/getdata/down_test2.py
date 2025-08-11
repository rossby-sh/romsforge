import os, random, time, yaml
from datetime import datetime as dt, timedelta
import cdsapi
import multiprocessing as mp

# 워커 (요청 + 리트라이)
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

# 하루 8타임
TIME8 = ["00:00","03:00","06:00","09:00","12:00","15:00","18:00","21:00"]

if __name__ == "__main__":
    # config 읽기
    with open("config.yml") as f:
        cfg = yaml.safe_load(f)

    start_dt = dt.strptime(cfg["time"]["start"], "%Y-%m-%d")  # "2025-06-30 12"
    end_dt   = dt.strptime(cfg["time"]["end"],   "%Y-%m-%d")  # "2025-08-02 00"
    pad_before = int(cfg["time"].get("pad_before_days", 0))
    pad_after  = int(cfg["time"].get("pad_after_days", 0))

    lat_min, lat_max = cfg["region"]["lat"]
    lon_min, lon_max = cfg["region"]["lon"]
    base_dir = cfg["output"]["base_dir"]

    tmp_dir = os.path.join(base_dir, "tmp_frc")
    os.makedirs(tmp_dir, exist_ok=True)

    req_start = (start_dt - timedelta(days=pad_before)).date()
    req_end   = (end_dt   + timedelta(days=pad_after)).date()

    s_date = start_dt.date()
    e_date = end_dt.date()

    # 시작 partial
    start_times = [h for h in TIME8 if int(h[:2]) >= start_dt.hour] or TIME8
    start_days  = [s_date]

    # 중간 full days
    middle_days = []
    # 패딩 앞쪽 full
    d = req_start
    while d < s_date:
        middle_days.append(d)
        d += timedelta(days=1)
    # 본구간 내부 full
    d = s_date + timedelta(days=1)
    while d < e_date:
        middle_days.append(d)
        d += timedelta(days=1)
    # 패딩 뒤쪽 full
    d = e_date + timedelta(days=1)
    while d <= req_end:
        middle_days.append(d)
        d += timedelta(days=1)
    # 시작/끝이 00시 또는 21시면 full로 처리
    if start_dt.hour == 0 and s_date not in start_days:
        middle_days.insert(0, s_date)
    if end_dt.hour == 21 and e_date not in middle_days:
        middle_days.append(e_date)

    # 끝 partial
    end_times = [h for h in TIME8 if int(h[:2]) <= end_dt.hour] or TIME8
    end_days  = [e_date]

    # 날짜 → 문자열 리스트 변환
    def ymd_lists(days_list):
        ys = sorted({d.strftime("%Y") for d in days_list})
        ms = sorted({d.strftime("%m") for d in days_list})
        ds = sorted({d.strftime("%d") for d in days_list})
        return ys, ms, ds

    dataset = "reanalysis-era5-single-levels"
    area = [lat_max, lon_min, lat_min, lon_max]

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

    jobs = []
    seg_info = []  # (kind, [seg_paths])
    period_tag = f"{req_start.strftime('%Y%m%d')}-{req_end.strftime('%Y%m%d')}"

    def add_segment(kind, vars_list, days_list, times_list, seg_name):
        if not days_list:
            return None
        ys, ms, ds = ymd_lists(days_list)
        req = {
            "product_type": ["reanalysis"],
            "variable": vars_list,
            "year": ys, "month": ms, "day": ds,
            "time": times_list,
            "format": "grib",
            "area": area,
        }
        seg_path = os.path.join(tmp_dir, f"{kind}_{seg_name}_{period_tag}.grib")
        if not os.path.exists(seg_path):
            jobs.append((dataset, req, seg_path))
        return seg_path

    for kind, vars_list in [("accum", accum_vars), ("inst", instant_vars)]:
        segs = []
        seg1 = add_segment(kind, vars_list, start_days, start_times, "seg1")
        if seg1: segs.append(seg1)
        seg2 = add_segment(kind, vars_list, middle_days, TIME8, "seg2")
        if seg2: segs.append(seg2)
        seg3 = add_segment(kind, vars_list, end_days, end_times, "seg3")
        if seg3: segs.append(seg3)
        seg_info.append((kind, segs))

    print(f"총 작업 {len(jobs)}개")
    if jobs:
        with mp.Pool(processes=6) as pool:
            pool.starmap(worker, jobs, chunksize=1)

    # 세그먼트 합쳐서 kind별 최종 파일 생성
#    for kind, segs in seg_info:
#        final_path = os.path.join(tmp_dir, f"{kind}_{period_tag}.grib")
#        if os.path.exists(final_path):
#            continue
#        with open(final_path, "wb") as w:
#            for seg in segs:
#                if os.path.exists(seg):
#                    with open(seg, "rb") as r:
#                        w.write(r.read())
#        print(f"[merge] {final_path}")
    # 세그먼트 합쳐서 kind별 최종 파일 생성 + 삭제
    for kind, segs in seg_info:
        final_path = os.path.join(tmp_dir, f"{kind}_{period_tag}.grib")

        # 머지 (최종 파일 없을 때만)
        if not os.path.exists(final_path):
            with open(final_path, "wb") as w:
                for seg in segs:
                    if os.path.exists(seg):
                        with open(seg, "rb") as r:
                            w.write(r.read())
            print(f"[merge] {final_path}")
        else:
            print(f"[skip merge] exists: {final_path}")

        # kind별 seg 삭제 (accum, inst 각각)
        removed = 0
        for seg in segs:
            if os.path.exists(seg):
                try:
                    os.remove(seg)
                    removed += 1
                except OSError as e:
                    print(f"[cleanup] 삭제 실패: {seg} ({e})")
        print(f"[cleanup] {kind}: removed {removed} seg files")


    print("--- Done ---")

