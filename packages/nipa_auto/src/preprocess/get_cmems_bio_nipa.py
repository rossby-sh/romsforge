# get_cmems_bio_nipa.py

from __future__ import annotations

import os
import json
from datetime import datetime, timedelta

import yaml
import numpy as np
import xarray as xr
import multiprocessing as mp
import copernicusmarine
from pathlib import Path

BASE = Path(__file__).resolve().parent
cfg_path = BASE / "config.yaml"


# -----------------------------------------------------------------------------
# read_cfg
# -----------------------------------------------------------------------------
def read_cfg(cfg_path):

    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f) or {}

    init_raw = cfg["initdate"]

    if isinstance(init_raw, datetime):
        init_dt = init_raw
    else:
        init_dt = datetime.fromisoformat(str(init_raw))

    init_dt = datetime(init_dt.year, init_dt.month, init_dt.day)

    cfg["time"] = {}
    cfg["time"]["init_dt"] = init_dt

    cfg.setdefault("workers", int(os.environ.get("CMEMS_WORKERS", "1")))

    cfg.setdefault("depth", {})
    cfg["depth"].setdefault("min", 0.0)
    cfg["depth"].setdefault("max", 5902.05810546875)

    base_dir = cfg["output"]["base_dir"].rstrip("/")

    cfg["output"]["raw_dir"]   = os.path.join(base_dir, "cmems_raw")
    cfg["output"]["cmems_dir"] = os.path.join(base_dir, "cmems")

    os.makedirs(cfg["output"]["raw_dir"], exist_ok=True)
    os.makedirs(cfg["output"]["cmems_dir"], exist_ok=True)

    cfg.setdefault("cmems_items", [

        {"name":"car",
         "dataset_id":"cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m",
         "variables":["dissic","ph","talk"]},

        {"name":"co2",
         "dataset_id":"cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m",
         "variables":["spco2"]},

        {"name":"nut",
         "dataset_id":"cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m",
         "variables":["fe","no3","po4","si"]},

        {"name":"pft",
         "dataset_id":"cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m",
         "variables":["chl","phyc"]},

        {"name":"bio",
         "dataset_id":"cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m",
         "variables":["o2","nppv"]},
    ])

    return cfg


# -----------------------------------------------------------------------------
# download_raw
# -----------------------------------------------------------------------------
def download_raw(cfg):

    init_dt = cfg["time"]["init_dt"]

    start_dt = init_dt - timedelta(days=1)
    end_dt   = init_dt + timedelta(days=1)

    start_iso = f"{start_dt:%Y-%m-%d}T00:00:00"
    end_iso   = f"{(end_dt + timedelta(days=1)):%Y-%m-%d}T00:00:00"

    tag = f"{start_dt:%Y%m%d}-{end_dt:%Y%m%d}"

    lon = tuple(cfg["region"]["lon"])
    lat = tuple(cfg["region"]["lat"])

    zmin = float(cfg["depth"]["min"])
    zmax = float(cfg["depth"]["max"])

    out_dir = cfg["output"]["raw_dir"]

    items = cfg["cmems_items"]


    def _download_one(item):

        out_name = f"CMEMS_{item['name']}_{tag}.nc"
        out_path = os.path.join(out_dir, out_name)

        if os.path.exists(out_path):
            return {"name":item["name"],"path":out_path,"status":"exists"}

        copernicusmarine.subset(

            dataset_id=item["dataset_id"],
            variables=item["variables"],

            minimum_longitude=lon[0],
            maximum_longitude=lon[1],

            minimum_latitude=lat[0],
            maximum_latitude=lat[1],

            start_datetime=start_iso,
            end_datetime=end_iso,

            minimum_depth=zmin,
            maximum_depth=zmax,

            output_filename=out_name,
            output_directory=out_dir,

            force_download=False,
        )

        return {"name":item["name"],"path":out_path,"status":"downloaded"}


    workers = int(cfg.get("workers",1))

    if workers <= 1:
        results = [_download_one(it) for it in items]
    else:
        with mp.Pool(processes=workers) as pool:
            results = list(pool.map(_download_one, items))


    meta = {
        "tag":tag,
        "start_iso":start_iso,
        "end_iso":end_iso,
        "created_at":datetime.utcnow().isoformat()+"Z",
        "results":results,
    }

    meta_path = os.path.join(out_dir,f"CMEMS_REQUEST_{tag}.json")

    with open(meta_path,"w") as f:
        json.dump(meta,f,indent=2)

    return results


# -----------------------------------------------------------------------------
# build_window_file
# -----------------------------------------------------------------------------
def build_window_file(cfg, raw_results):

    init_dt = cfg["time"]["init_dt"]

    start_dt = init_dt - timedelta(days=1)
    end_dt   = init_dt + timedelta(days=1)

    t0 = np.datetime64(start_dt.date())
    t1 = np.datetime64(end_dt.date())

    parts = []

    for r in raw_results:

        ds = xr.open_dataset(r["path"],decode_times=True)

        if "time" not in ds.coords:

            for nm,da in ds.coords.items():

                if np.issubdtype(da.dtype,np.datetime64):
                    ds = ds.rename({nm:"time"})
                    break


        ds = ds.resample(time="1D").mean()

        ds = ds.sel(time=slice(t0,t1))

        parts.append(ds)


    merged = xr.merge(parts,compat="override")

    merged = merged.sortby("time")

    out_file = os.path.join(
        cfg["output"]["cmems_dir"],
        f"cmems_bio_{init_dt:%Y%m%d}.nc"
    )

    tmp = out_file + ".tmp"

    merged.to_netcdf(tmp)

    os.replace(tmp,out_file)

    print("[OK] created",out_file)


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():

    cfg = read_cfg(cfg_path)

    raw_results = download_raw(cfg)

    build_window_file(cfg, raw_results)

    print("[OK] raw_dir  =",cfg["output"]["raw_dir"])
    print("[OK] cmems_dir=",cfg["output"]["cmems_dir"])


if __name__ == "__main__":
    main()
