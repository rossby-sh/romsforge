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
    # config 읽기 (날짜만 줘도 됨: %Y-%m-%d)
    with open("config.yml") as f:
        cfg = yaml.safe_load(f)

    # 날짜만 입력해도 동작하게: 시는 00시로 들어감
    start_dt = dt.strptime(cfg["time"]["start"], "%Y-%m-%d")
    end_dt   = dt.strptime(cfg["time"]["end"],   "%Y-%m-%d")
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

    # 시작 partial/끝 partial → 날짜만 줘도 항상 생성되도록
    start_times = [h for h in TIME8 if int(h[:2]) >= start_dt.hour] or TIME8
    start_days  = [s_date]

    end_times = [h for h in TIME8 if int(h[:2]) <= end_dt.hour] or TIME8
    end_days  = [e_date]

    # 중간 full days(패딩 포함)
    middle_days = []
    d = req_start
    while d < s_date:
        middle_days.append(d); d += timedelta(days=1)
    d = s_date + timedelta(days=1)
    while d < e_date:
        middle_days.append(d); d += timedelta(days=1)
    d = e_date + timedelta(days=1)
    while d <= req_end:
        middle_days.append(d); d += timedelta(days=1)

    # 날짜 → 문자열 리스트 변환
    def ymd_lists(days_list):
        ys = sorted({d.strftime("%Y") for d in days_list})
        ms = sorted({d.strftime("%m") for d in days_list})
        ds = sorted({d.strftime("%d") for d in days_list})
        return ys, ms, ds

    # 데이터셋
    dataset_sl = "reanalysis-era5-single-levels"
    dataset_pl = "reanalysis-era5-pressure-levels"  # ← 1000hPa specific humidity
    area = [lat_max, lon_min, lat_min, lon_max]

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
    # pressure-level: 1000hPa specific humidity
    q1000_vars = ["specific_humidity"]
    q1000_extra = {"pressure_level": ["1000"]}

    jobs = []
    # kind별 세그먼트 경로 저장해서 나중에 병합+삭제
    seg_info = []  # (kind, seg_paths)
    period_tag = f"{req_start.strftime('%Y%m%d')}-{req_end.strftime('%Y%m%d')}"

    def add_segment(dataset_name, kind, vars_list, days_list, times_list, seg_name, extra=None):
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
        if extra:
            req.update(extra)  # pressure_level 등 추가
        seg_path = os.path.join(tmp_dir, f"{kind}_{seg_name}_{period_tag}.grib")
        if not os.path.exists(seg_path):
            jobs.append((dataset_name, req, seg_path))
        return seg_path

    # 1) accum (single-levels)
    segs_accum = []
    s = add_segment(dataset_sl, "accum", accum_vars, start_days,  start_times, "seg1")
    if s: segs_accum.append(s)
    s = add_segment(dataset_sl, "accum", accum_vars, middle_days, TIME8,       "seg2")
    if s: segs_accum.append(s)
    s = add_segment(dataset_sl, "accum", accum_vars, end_days,    end_times,   "seg3")
    if s: segs_accum.append(s)
    seg_info.append(("accum", segs_accum))

    # 2) inst (single-levels)
    segs_inst = []
    s = add_segment(dataset_sl, "inst", instant_vars, start_days,  start_times, "seg1")
    if s: segs_inst.append(s)
    s = add_segment(dataset_sl, "inst", instant_vars, middle_days, TIME8,       "seg2")
    if s: segs_inst.append(s)
    s = add_segment(dataset_sl, "inst", instant_vars, end_days,    end_times,   "seg3")
    if s: segs_inst.append(s)
    seg_info.append(("inst", segs_inst))

    # 3) q1000 (pressure-levels, specific humidity @1000hPa)
    segs_q1000 = []
    s = add_segment(dataset_pl, "q1000", q1000_vars, start_days,  start_times, "seg1", extra=q1000_extra)
    if s: segs_q1000.append(s)
    s = add_segment(dataset_pl, "q1000", q1000_vars, middle_days, TIME8,       "seg2", extra=q1000_extra)
    if s: segs_q1000.append(s)
    s = add_segment(dataset_pl, "q1000", q1000_vars, end_days,    end_times,   "seg3", extra=q1000_extra)
    if s: segs_q1000.append(s)
    seg_info.append(("q1000", segs_q1000))

    print(f"총 작업 {len(jobs)}개")
    if jobs:
        with mp.Pool(processes=9) as pool:   # 원하면 2~3으로 낮춰도 됨
            pool.starmap(worker, jobs, chunksize=1)

    # 머지 + seg 삭제 (accum, inst, q1000 각각)
    for kind, segs in seg_info:
        final_path = os.path.join(tmp_dir, f"{kind}_{period_tag}.grib")
        # merge (없을 때만 생성)
        if not os.path.exists(final_path):
            with open(final_path, "wb") as w:
                for seg in segs:
                    if os.path.exists(seg):
                        with open(seg, "rb") as r:
                            w.write(r.read())
            print(f"[merge] {final_path}")
        else:
            print(f"[skip merge] exists: {final_path}")
        # cleanup
        removed = 0
        for seg in segs:
            if os.path.exists(seg):
                try:
                    os.remove(seg); removed += 1
                except OSError as e:
                    print(f"[cleanup] 삭제 실패: {seg} ({e})")
        print(f"[cleanup] {kind}: removed {removed} seg files")

    print("--- Done ---")

