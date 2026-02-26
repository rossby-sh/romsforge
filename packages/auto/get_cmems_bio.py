
# cmems_bio_pipeline_fortran_style.py
from __future__ import annotations

import os
import json
from datetime import datetime, timedelta

import yaml
import numpy as np
import xarray as xr
import multiprocessing as mp
import copernicusmarine


# -----------------------------------------------------------------------------
# [SUB] read_cfg
# -----------------------------------------------------------------------------
def read_cfg(cfg_path: str) -> dict:
    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f) or {}

    # required
    cfg["time"]["start_dt"] = datetime.strptime(cfg["time"]["start"], "%Y-%m-%d")
    cfg["time"]["end_dt"]   = datetime.strptime(cfg["time"]["end"],   "%Y-%m-%d")

    # defaults
    cfg.setdefault("workers", int(os.environ.get("CMEMS_WORKERS", "1")))
    cfg.setdefault("daily_collapse", "mean")     # mean | first
    cfg.setdefault("force_shard", False)

    cfg.setdefault("depth", {})
    cfg["depth"].setdefault("min", 0.0)
    cfg["depth"].setdefault("max", 5902.05810546875)

    cfg.setdefault("output", {})
    base_dir = cfg["output"]["base_dir"].rstrip("/")
    cfg["output"]["raw_dir"]   = os.path.join(base_dir, "cmems_bio_raw")
    cfg["output"]["daily_dir"] = os.path.join(base_dir, "cmems_bio_daily")
    os.makedirs(cfg["output"]["raw_dir"], exist_ok=True)
    os.makedirs(cfg["output"]["daily_dir"], exist_ok=True)

    # items (allow config)
    cfg.setdefault("cmems_items", [
        {"name": "car", "dataset_id": "cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m", "variables": ["dissic", "ph", "talk"]},
        {"name": "co2", "dataset_id": "cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m", "variables": ["spco2"]},
        {"name": "nut", "dataset_id": "cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m", "variables": ["fe", "no3", "po4", "si"]},
        {"name": "pft", "dataset_id": "cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m", "variables": ["chl", "phyc"]},
        {"name": "bio", "dataset_id": "cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m", "variables": ["o2", "nppv"]},
    ])

    # optional: rename map (cmems var -> your ogcm std var)
    # You said you will add:
    # bio_ogcm_var_name: { NO3: no3, PO4: po4, oxygen: o2, chl: chl, phyt: phyc, ... }
    cfg.setdefault("bio_ogcm_var_name", {})

    return cfg


# -----------------------------------------------------------------------------
# [SUB] download_raw
# -----------------------------------------------------------------------------
def download_raw(cfg: dict) -> list[dict]:
    start_dt: datetime = cfg["time"]["start_dt"]
    end_dt: datetime   = cfg["time"]["end_dt"]

    # request window: end + 1 day (00:00) to include end date safely
    start_iso = f"{start_dt:%Y-%m-%d}T00:00:00"
    end_iso   = f"{(end_dt + timedelta(days=1)):%Y-%m-%d}T00:00:00"
    tag = f"{start_dt:%Y%m%d}-{end_dt:%Y%m%d}"

    lon = tuple(cfg["region"]["lon"])
    lat = tuple(cfg["region"]["lat"])
    zmin = float(cfg["depth"]["min"])
    zmax = float(cfg["depth"]["max"])

    out_dir = cfg["output"]["raw_dir"]
    items = cfg["cmems_items"]

    def _download_one(item: dict) -> dict:
        out_name = f"CMEMS_{item['name']}_{tag}.nc"
        out_path = os.path.join(out_dir, out_name)
        if os.path.exists(out_path):
            return {"name": item["name"], "path": out_path, "status": "exists"}

        copernicusmarine.subset(
            dataset_id=item["dataset_id"],
            variables=item["variables"],
            minimum_longitude=lon[0], maximum_longitude=lon[1],
            minimum_latitude=lat[0],  maximum_latitude=lat[1],
            start_datetime=start_iso, end_datetime=end_iso,
            minimum_depth=zmin, maximum_depth=zmax,
            output_filename=out_name,
            output_directory=out_dir,
            force_download=False,
        )
        return {"name": item["name"], "path": out_path, "status": "downloaded"}

    workers = int(cfg.get("workers", 1))
    if workers <= 1:
        results = [_download_one(it) for it in items]
    else:
        with mp.Pool(processes=workers) as pool:
            results = list(pool.map(_download_one, items))

    # write request meta
    meta = {
        "tag": tag,
        "start_iso": start_iso,
        "end_iso": end_iso,
        "lon": list(lon),
        "lat": list(lat),
        "depth_min": zmin,
        "depth_max": zmax,
        "items": items,
        "results": results,
        "created_at": datetime.utcnow().isoformat() + "Z",
    }
    meta_path = os.path.join(out_dir, f"CMEMS_REQUEST_{tag}.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    cfg["time"]["tag"] = tag
    return results


# -----------------------------------------------------------------------------
# [SUB] make_daily_shards
# -----------------------------------------------------------------------------
def make_daily_shards(cfg: dict, raw_results: list[dict]) -> None:
    start_dt: datetime = cfg["time"]["start_dt"]
    end_dt: datetime   = cfg["time"]["end_dt"]
    t0 = np.datetime64(start_dt.date())
    t1 = np.datetime64(end_dt.date())

    out_dir = cfg["output"]["daily_dir"]
    tag = cfg["time"]["tag"]

    collapse = str(cfg.get("daily_collapse", "mean")).lower()  # mean|first
    force = bool(cfg.get("force_shard", False))

    # rename map: cmems var -> std var (we build from bio_ogcm_var_name inverse)
    # You provide std->cmems mapping, e.g. NO3:no3. We need cmems->std for rename.

    def _open_and_daily(path: str) -> xr.Dataset:
        ds = xr.open_dataset(path, decode_times=True, mask_and_scale=True)
        # time coord name unify
        if "time" not in ds.coords:
            # best-effort: find CF axis T
            tname = None
            for nm, da in ds.coords.items():
                if str(da.attrs.get("axis", "")).upper() == "T":
                    tname = nm
                    break
            if tname is None:
                # last resort: any datetime coord
                for nm, da in ds.coords.items():
                    if np.issubdtype(da.dtype, np.datetime64):
                        tname = nm
                        break
            if tname is None:
                raise RuntimeError(f"cannot find time coord in {path}")
            ds = ds.rename({tname: "time"})

        # daily collapse
        day = ds["time"].astype("datetime64[D]").astype("datetime64[ns]")
        ds = ds.assign_coords(day=("time", day.data))
        if collapse == "mean":
            ds = ds.groupby("day").mean("time", keep_attrs=True)
        elif collapse == "first":
            ds = ds.groupby("day").first("time", keep_attrs=True)
        else:
            raise ValueError(f"unknown daily_collapse: {collapse}")
        ds = ds.rename({"day": "time"})

        # rename variables to std names (cmems -> std)
        # only rename those present

        return ds

    daily_sets = []
    for r in raw_results:
        p = r["path"]
        if not os.path.exists(p):
            raise FileNotFoundError(p)
        daily_sets.append(_open_and_daily(p))

    # union time list
    all_times = None
    for ds in daily_sets:
        tt = ds["time"].values
        all_times = tt if all_times is None else np.union1d(all_times, tt)

    # filter to [t0, t1]
    days = []
    for tt in all_times:
        d = tt.astype("datetime64[D]")
        if d < t0 or d > t1:
            continue
        days.append(d.astype("datetime64[ns]"))
    days = np.array(sorted(set(days)), dtype="datetime64[ns]")

    # write each day merged
    for day in days:
        out_path = os.path.join(out_dir, f"cmems_bio_{str(day)[:10].replace('-','')}.nc")
        if os.path.exists(out_path) and (not force):
            continue

        parts = []
        for ds in daily_sets:
            if day in ds["time"].values:
                parts.append(ds.sel(time=[day]))

        if not parts:
            continue

        merged = xr.merge(parts, compat="override", combine_attrs="override")

        # ensure time is a dimension with length=1
        if ("time" in merged.coords) and ("time" not in merged.dims):
            merged = merged.expand_dims("time")

        merged.attrs["source"] = "CMEMS subset -> daily merged shard"
        merged.attrs["day"] = str(day)[:10]
        merged.attrs["request_tag"] = tag


        tmp = out_path + ".tmp"
        merged.to_netcdf(tmp)
        os.replace(tmp, out_path)


# -----------------------------------------------------------------------------
# [MAIN]
# -----------------------------------------------------------------------------
def main():
    cfg = read_cfg("config.yaml")

    # (1) download raw
    raw_results = download_raw(cfg)

    # (2) postprocess to daily shards (merged 1 file per day)
    make_daily_shards(cfg, raw_results)

    print(f"[OK] raw_dir  = {cfg['output']['raw_dir']}")
    print(f"[OK] daily_dir= {cfg['output']['daily_dir']}")


if __name__ == "__main__":
    main()
