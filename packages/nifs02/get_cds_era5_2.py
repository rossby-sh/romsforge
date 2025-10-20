
# inspect_era5_time.py
import xarray as xr
import numpy as np

def _summ_coord(ds, name):
    if name not in ds.coords:
        return f"{name}: (없음)"
    da = ds.coords[name]
    vals = da.values
    # 스칼라(0-d)
    if not hasattr(vals, "ndim") or vals.ndim == 0:
        return f"{name}: {vals} (scalar)"
    # 1-d 이상
    return f"{name}: {vals[0]}  →  {vals[-1]}  (len={vals.shape[0]})"

def _summ_lonlat(ds):
    # cfgrib은 보통 latitude/longitude 이름 사용
    lon_name = "longitude" if "longitude" in ds.coords else ("lon" if "lon" in ds.coords else None)
    lat_name = "latitude"  if "latitude"  in ds.coords else ("lat" if "lat" in ds.coords else None)
    out = []
    if lon_name is not None:
        v = ds.coords[lon_name].values
        if v.ndim == 1:
            out.append(f"lon: {float(v[0])} → {float(v[-1])} (n={v.size})")
        else:
            out.append(f"lon: shape={v.shape}")
    else:
        out.append("lon: (없음)")
    if lat_name is not None:
        v = ds.coords[lat_name].values
        if v.ndim == 1:
            out.append(f"lat: {float(v[0])} → {float(v[-1])} (n={v.size})")
        else:
            out.append(f"lat: shape={v.shape}")
    else:
        out.append("lat: (없음)")
    return " | ".join(out)

def inspect_grib(path):
    print(f"\n=== {path} ===")
    try:
        ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
        print(_summ_coord(ds, "time"))
        print(_summ_coord(ds, "valid_time"))
        print(_summ_coord(ds, "step"))
        print(_summ_lonlat(ds))
        print("coords:", list(ds.coords))
        ds.close()
    except Exception as e:
        print("단일 dataset으로 못 열림:", e)
        print("→ cfgrib.open_datasets로 서브-데이터셋들 확인")
        from cfgrib import open_datasets
        dsets = open_datasets(path, indexpath="")
        for i, sub in enumerate(dsets):
            print(f"\n-- Sub-dataset {i} --")
            print(_summ_coord(sub, "time"))
            print(_summ_coord(sub, "valid_time"))
            print(_summ_coord(sub, "step"))
            print(_summ_lonlat(sub))
            print("coords:", list(sub.coords))

if __name__ == "__main__":
    accum_path = "/home/shjo/data/nifs02/sep/era5/accum_20250831-20251004.grib"
    inst_path  = "/home/shjo/data/nifs02/sep/era5/inst_20250831-20251004.grib"
    inspect_grib(inst_path)
    print("!!!!!!!!!!!!!!!!")
    inspect_grib(accum_path)
