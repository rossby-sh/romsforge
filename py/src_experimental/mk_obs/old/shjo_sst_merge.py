import xarray as xr
import numpy as np
from netCDF4 import Dataset

pth='D:/shjo/ROMS_inputs/obs/'

# 입력 파일 경로 (수정해도 됨)
sst_file = pth+"obs_SST_OSTIA_30km.nc"
chl_file = pth+"obs_phyt_27km.nc" 
merged_file = pth+"obs_sst_chl_30km.nc"

# 파일 열기
ds_sst = xr.open_dataset(sst_file,decode_times=False)
ds_chl = xr.open_dataset(chl_file,decode_times=False)
ds_sst["obs_time"].values = np.floor(ds_sst["obs_time"].values + 1e-6)

ds_chl["obs_time"].values

if "survey" in ds_sst.dims:
    ds_sst = ds_sst.isel(survey=0).drop_vars("survey_time", errors="ignore")
if "survey" in ds_chl.dims:
    ds_chl = ds_chl.isel(survey=0).drop_vars("survey_time", errors="ignore")

# obs_type 명시적으로 재설정
ds_sst["obs_type"] = xr.DataArray(np.full(ds_sst.dims["datum"], 6, dtype=np.int32), dims="datum")
ds_chl["obs_type"] = xr.DataArray(np.full(ds_chl.dims["datum"], 10, dtype=np.int32), dims="datum")

# 병합
merged = xr.concat([ds_sst, ds_chl], dim="datum")

# obs_time 기준 정렬
sorted_idx = np.argsort(merged["obs_time"].values)

merged = merged.isel(datum=sorted_idx)

# survey_time, Nobs 계산
obs_time = merged["obs_time"].values

survey_time, inverse = np.unique(obs_time, return_inverse=True)
Nobs = np.bincount(inverse)

# survey 차원 추가 (빈 차원이 아닌, 메타데이터용)
merged["survey_time"] = xr.DataArray(survey_time, dims="survey")
merged["Nobs"] = xr.DataArray(Nobs, dims="survey")

# obs_variance 보존
if "obs_variance" in ds_sst:
    merged["obs_variance"] = ds_sst["obs_variance"]




# merged.to_netcdf(merged_file)



with Dataset(merged_file, "w", format="NETCDF4") as nc:
    # Dimensions
    nc.createDimension("datum", merged.dims["datum"])
    nc.createDimension("survey", merged.dims["survey"])
    nc.createDimension("state_variable", merged.dims["state_variable"])

    # datum variables
    for var in ["obs_time", "obs_value", "obs_error", "obs_type", 
                "obs_lat", "obs_lon", "obs_depth", "obs_Xgrid", "obs_Ygrid", "obs_Zgrid"]:
        if var in merged:
            dtype = "f4" if merged[var].dtype.kind == "f" else "i4"
            v = nc.createVariable(var, dtype, ("datum",))
            v[:] = merged[var].values.astype(dtype)

    # survey variables
    nc.createVariable("survey_time", "f4", ("survey",))[:] = merged["survey_time"].values.astype("f4")
    nc.createVariable("Nobs", "i4", ("survey",))[:] = merged["Nobs"].values.astype("i4")

    # state_variable
    nc.createVariable("obs_variance", "f4", ("state_variable",))[:] = merged["obs_variance"].values.astype("f4")


print(f"✅ 병합 완료: {merged_file}")



