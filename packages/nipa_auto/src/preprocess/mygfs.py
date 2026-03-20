#!/usr/bin/env python3
from __future__ import annotations

import os
import yaml
import requests
import subprocess
import numpy as np
import xarray as xr
from datetime import datetime as dt, timedelta
from pathlib import Path

BASE = Path(__file__).resolve().parent
cfg_path = BASE / "config.yaml"

S3_BASE = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
RES = "0p25"

RUN_HOURS = (0, 6, 12, 18)

TIME_REF = dt(2000,1,1)
STEP_HOURS = 3


# =========================
# helpers
# =========================

def as_dt(x):
    if isinstance(x, dt):
        return x
    return dt.strptime(x,"%Y-%m-%d %H:%M:%S")


def make_regular_times(t_st, t_ed, step):

    step = timedelta(hours=step)
    cur = t_st

    out = []

    while cur <= t_ed:
        out.append(cur)
        cur += step

    return out


def floor_to_run(t):

    cands = [
        t.replace(hour=h,minute=0,second=0,microsecond=0)
        for h in RUN_HOURS
    ]

    past = [x for x in cands if x <= t]

    if past:
        return max(past)

    return (t - timedelta(days=1)).replace(hour=18)


def fhr_from(run, valid):

    return int((valid - run).total_seconds() / 3600)


def gfs_url(run, fhr):

    ymd = run.strftime("%Y%m%d")
    hh = run.strftime("%H")

    return f"{S3_BASE}/gfs.{ymd}/{hh}/atmos/gfs.t{hh}z.pgrb2.{RES}.f{fhr:03d}"


def download(url, out):

    if os.path.exists(out):
        return

    print("download:", url)

    r = requests.get(url, stream=True)

    r.raise_for_status()

    with open(out,"wb") as f:
        for c in r.iter_content(1024*1024):
            f.write(c)


def run(cmd, inp=None):

    p = subprocess.run(
        cmd,
        input=inp,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if p.returncode != 0:
        raise RuntimeError(p.stderr)

    return p.stdout


def detect_latlon_names(ds):

    for cand in ("lat","latitude","Latitude"):
        if cand in ds.coords:
            lat_name = cand
            break
    else:
        raise RuntimeError("latitude coord not found")

    for cand in ("lon","longitude","Longitude"):
        if cand in ds.coords:
            lon_name = cand
            break
    else:
        raise RuntimeError("longitude coord not found")

    return lat_name, lon_name


# =========================
# main
# =========================

def main():

    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t_st = as_dt(cfg["bry_start_date"])
    t_ed = as_dt(cfg["bry_end_date"])

    lat_rng = tuple(cfg["region"]["lat"])
    lon_rng = tuple(cfg["region"]["lon"])

    base_dir = cfg["output"]["base_dir"]

    gfs_dir = os.path.join(base_dir,"gfs")
    raw_dir = os.path.join(gfs_dir,"raw")

    os.makedirs(raw_dir, exist_ok=True)

    print("[INFO] window:", t_st, "→", t_ed)

    target_times = make_regular_times(t_st, t_ed, STEP_HOURS)

    pieces = []

    for t in target_times:

        run_inst = floor_to_run(t)
        fhr_inst = fhr_from(run_inst,t)

        t_prev = t - timedelta(hours=3)

        run_cum = floor_to_run(t_prev)
        fhr_cum = fhr_from(run_cum,t)

        print("\nTARGET",t)
        print("instant:",run_inst,fhr_inst)
        print("cumulative:",run_cum,fhr_cum)

        url_i = gfs_url(run_inst,fhr_inst)
        url_c = gfs_url(run_cum,fhr_cum)

        raw_i = os.path.join(
            raw_dir,
            f"{run_inst:%Y%m%d%H}_f{fhr_inst:03d}.grib2"
        )

        raw_c = os.path.join(
            raw_dir,
            f"{run_cum:%Y%m%d%H}_f{fhr_cum:03d}.grib2"
        )

        download(url_i,raw_i)
        download(url_c,raw_c)

        thin = os.path.join(raw_dir,"thin.grib2")

        if os.path.exists(thin):
            os.remove(thin)

        # instant vars
        inv_i = run(["wgrib2",raw_i]).splitlines()

        want_i = [
        ":UGRD:10 m above ground:",
        ":VGRD:10 m above ground:",
        ":TMP:2 m above ground:",
        ":SPFH:2 m above ground:",
        ":PRMSL:mean sea level:"
        ]

        sel = []

        for pat in want_i:
            for ln in inv_i:
                if pat in ln:
                    sel.append(ln)
                    break

        run(
        ["wgrib2",raw_i,"-i","-grib_out",thin],
        "\n".join(sel)+"\n"
        )

        # cumulative vars
        inv_c = run(["wgrib2",raw_c]).splitlines()

        want_c = [
        ":DSWRF:surface:",
        ":DLWRF:surface:",
        ":APCP:surface:"
        ]

        sel = []

        for pat in want_c:
            for ln in inv_c:
                if pat in ln:
                    sel.append(ln)
                    break

        run(
        ["wgrib2",raw_c,"-i","-append","-grib_out",thin],
        "\n".join(sel)+"\n"
        )

        tmp_nc = os.path.join(raw_dir,"tmp.nc")

        run(["wgrib2",thin,"-netcdf",tmp_nc])

        ds = xr.open_dataset(tmp_nc)

        lat_name, lon_name = detect_latlon_names(ds)

        ds = ds.sel({
            lat_name: slice(*lat_rng),
            lon_name: slice(*lon_rng)
        })

        days = (t - TIME_REF).total_seconds() / 86400

        if ds.sizes.get("time",0) > 1:
            ds = ds.isel(time=0)

        ds = ds.assign_coords(time=("time",[days]))

        ds["time"].attrs["units"] = "days since 2000-01-01 00:00:00"
        ds["time"].attrs["calendar"] = "proleptic_gregorian"

        pieces.append(ds)

    final = xr.concat(pieces, dim="time")

    out = os.path.join(gfs_dir,"gfs_korea.nc")

    final.to_netcdf(out)

    print("\nsaved:",out)


if __name__ == "__main__":
    main()
