
#!/usr/bin/env python3
import os
import re
import glob
import datetime as dt

import numpy as np
import xarray as xr

# --- config ------------------------------------------------------------
SRC_GLOB = "/home/shjo/data/nifs02/sep_isl_test/modis/AQUA_MODIS.*.NW_Pacific.chlor_a_NRT.nc"
PATTERN = re.compile(r"AQUA_MODIS\.(\d{8})\.NW_Pacific\.chlor_a_NRT\.nc$")

TIME_UNITS = "days since 2000-01-01 00:00:00"
REF = dt.datetime(2000, 1, 1)

# 원본 파일을 직접 덮어씀
OVERWRITE = True
# ---------------------------------------------------------------------


def yyyymmdd_to_days(yyyymmdd: str) -> float:
    d = dt.datetime.strptime(yyyymmdd, "%Y%m%d")
    return (d - REF).total_seconds() / 86400.0


def fix_one(path: str, overwrite: bool = True) -> str:
    fname = os.path.basename(path)
    m = PATTERN.match(fname)
    if not m:
        raise ValueError(f"Filename does not match pattern: {fname}")

    yyyymmdd = m.group(1)
    tval = yyyymmdd_to_days(yyyymmdd)

    # open (decode_times 끔: 우리가 직접 time 관리)
    ds = xr.open_dataset(path, decode_times=False)

    # sanity checks
    if "time" in ds.dims or "time" in ds.variables:
        ds.close()
        raise RuntimeError(f"'time' already exists in: {fname}")

    if "chlor_a" not in ds.variables:
        ds.close()
        raise RuntimeError(f"'chlor_a' not found in: {fname}")

    if ds["chlor_a"].ndim != 2:
        ds.close()
        raise RuntimeError(
            f"Expected chlor_a(lat, lon) 2D, got ndim={ds['chlor_a'].ndim} in: {fname}"
        )

    # add time dimension
    ds2 = ds.expand_dims(time=[np.float64(tval)])
    ds2["time"].attrs["units"] = TIME_UNITS
    ds2["time"].attrs["calendar"] = "standard"

    # overwrite with atomic replace
    tmp_path = path + ".tmp"
    ds2.to_netcdf(tmp_path)

    ds.close()
    ds2.close()

    os.replace(tmp_path, path)
    return path


def main():
    files = sorted(glob.glob(SRC_GLOB))
    if not files:
        raise SystemExit(f"No files matched: {SRC_GLOB}")

    for f in files:
        try:
            out = fix_one(f, overwrite=OVERWRITE)
            print(f"[OK] overwritten: {out}")
        except Exception as e:
            print(f"[FAIL] {f}: {e}")


if __name__ == "__main__":
    main()
