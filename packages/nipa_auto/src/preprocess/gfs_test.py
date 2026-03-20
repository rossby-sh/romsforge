#!/usr/bin/env python3
from __future__ import annotations

import os
import yaml
import numpy as np
import xarray as xr
from datetime import datetime as dt, timedelta
from typing import Any, Dict
from pathlib import Path

BASE = Path(__file__).resolve().parent
cfg_path = BASE / "config.yaml"

TIME_REF = dt(2000,1,1)

RUN_HOURS = (0,6,12,18)

# -----------------------
# helpers
# -----------------------

def as_dt(x: Any) -> dt:
    if isinstance(x, dt):
        return x
    return dt.strptime(x,"%Y-%m-%d %H:%M:%S")


def make_regular_times(t_st:dt, t_ed:dt, step:int):
    cur=t_st
    step=timedelta(hours=step)
    out=[]
    while cur<=t_ed:
        out.append(cur)
        cur+=step
    return out


def floor_to_run(t:dt)->dt:
    cands=[t.replace(hour=h,minute=0,second=0,microsecond=0) for h in RUN_HOURS]
    past=[x for x in cands if x<=t]
    if past:
        return max(past)
    return (t-timedelta(days=1)).replace(hour=18,minute=0,second=0,microsecond=0)


def fhr_from(run:dt, valid:dt)->int:
    return int((valid-run).total_seconds()/3600)


# -----------------------
# variable mapping
# -----------------------

VARMAP={
"u10":"UGRD_10maboveground",
"v10":"VGRD_10maboveground",
"tair":"TMP_2maboveground",
"qair":"SPFH_2maboveground",
"pair":"PRMSL_meansealevel",
"swrad":"DSWRF_surface",
"lwrad":"DLWRF_surface",
"rain":"APCP_surface"
}


# -----------------------
# main
# -----------------------

def main():

    with open(cfg_path) as f:
        cfg=yaml.safe_load(f)

    t_st=as_dt(cfg["bry_start_date"])
    t_ed=as_dt(cfg["bry_end_date"])

    lat_rng=tuple(cfg["region"]["lat"])
    lon_rng=tuple(cfg["region"]["lon"])

    base_dir=cfg["output"]["base_dir"]
    gfs_dir=os.path.join(base_dir,"gfs")
    os.makedirs(gfs_dir,exist_ok=True)

    target_times=make_regular_times(t_st,t_ed,3)

    pieces=[]

    for t in target_times:

        run_inst=floor_to_run(t)
        fhr_inst=fhr_from(run_inst,t)

        t_prev=t-timedelta(hours=3)
        run_cum=floor_to_run(t_prev)
        fhr_cum=fhr_from(run_cum,t)

        print("target:",t)
        print(" instant:",run_inst,fhr_inst)
        print(" cumulative:",run_cum,fhr_cum)

        # ----------------
        # open GFS files
        # ----------------

        fn_inst=f"gfs_{run_inst:%Y%m%d%H}_f{fhr_inst:03d}.nc"
        fn_cum =f"gfs_{run_cum:%Y%m%d%H}_f{fhr_cum:03d}.nc"

        ds_i=xr.open_dataset(os.path.join(gfs_dir,fn_inst))
        ds_c=xr.open_dataset(os.path.join(gfs_dir,fn_cum))

        # ----------------
        # subset region
        # ----------------

        ds_i=ds_i.sel(lat=slice(*lat_rng),lon=slice(*lon_rng))
        ds_c=ds_c.sel(lat=slice(*lat_rng),lon=slice(*lon_rng))

        # ----------------
        # variable build
        # ----------------

        out=xr.Dataset()

        out["u10"]=ds_i[VARMAP["u10"]]
        out["v10"]=ds_i[VARMAP["v10"]]
        out["tair"]=ds_i[VARMAP["tair"]]
        out["qair"]=ds_i[VARMAP["qair"]]
        out["pair"]=ds_i[VARMAP["pair"]]

        out["swrad"]=ds_c[VARMAP["swrad"]]
        out["lwrad"]=ds_c[VARMAP["lwrad"]]
        out["rain"]=ds_c[VARMAP["rain"]]

        # ----------------
        # time coord
        # ----------------

        days=(t-TIME_REF).total_seconds()/86400

        out=out.expand_dims(time=[days])
        out.time.attrs["units"]="days since 2000-01-01 00:00:00"
        out.time.attrs["calendar"]="proleptic_gregorian"

        pieces.append(out)

    # -----------------------
    # concat
    # -----------------------

    final=xr.concat(pieces,dim="time")

    out_path=os.path.join(gfs_dir,"gfs_forcing.nc")

    final.to_netcdf(out_path)

    print("saved:",out_path)


if __name__=="__main__":
    main()

