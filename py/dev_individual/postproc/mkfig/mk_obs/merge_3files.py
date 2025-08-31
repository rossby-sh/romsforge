import xarray as xr
import numpy as np
from netCDF4 import Dataset

pth = 'D:/shjo/ROMS_inputs/obs/pro/'

# 입력 파일들
sst_file = pth + "obs_SST_OSTIA_30km.nc"
chl_file = pth + "obs_phyt_27km.nc"
kodc_file = pth + "obs_KODC.nc"
merged_file = pth + "test_merge.nc"

# 열기
ds_sst = xr.open_dataset(sst_file, decode_times=False)
ds_chl = xr.open_dataset(chl_file, decode_times=False)
ds_kodc = xr.open_dataset(kodc_file, decode_times=False)

# 시간 정리 (float64 rounding 방지)
for ds in [ds_sst, ds_chl]:
    ds["obs_time"].values[:] = np.floor(ds["obs_time"].values + 1e-6)

# survey 차원 제거 함수
def drop_survey(ds):
    if "survey" in ds.dims:
        ds = ds.isel(survey=0).drop_vars("survey_time", errors="ignore")
        ds = ds.drop_dims("survey", errors="ignore")
    return ds

ds_sst = drop_survey(ds_sst)
ds_chl = drop_survey(ds_chl)
ds_kodc = drop_survey(ds_kodc)

# obs_type 재설정 (sst=6, chl=9), kodc는 그대로
ds_sst["obs_type"] = xr.DataArray(np.full(ds_sst.dims["datum"], 6, dtype=np.int32), dims="datum")
ds_chl["obs_type"] = xr.DataArray(np.full(ds_chl.dims["datum"], 11, dtype=np.int32), dims="datum")

ds_sst["obs_provenance"] = xr.DataArray(np.full(ds_sst.dims["datum"], 2, dtype=np.int32), dims="datum")
ds_chl["obs_provenance"] = xr.DataArray(np.full(ds_chl.dims["datum"], 14, dtype=np.int32), dims="datum")

# 병합
merged = xr.concat([ds_sst, ds_chl, ds_kodc], dim="datum")

# 정렬
sorted_idx = np.argsort(merged["obs_time"].values)
merged = merged.isel(datum=sorted_idx)

# survey 계산
obs_time = merged["obs_time"].values
survey_time, inverse = np.unique(obs_time, return_inverse=True)
Nobs = np.bincount(inverse)
merged["survey_time"] = xr.DataArray(survey_time, dims="survey")
merged["Nobs"] = xr.DataArray(Nobs, dims="survey")

# obs_variance 계산 (obs_type별 평균 error^2)
# obs_variance = np.ones(12, dtype="f4")
# for var_type in np.unique(merged["obs_type"].values):
#     idx = merged["obs_type"].values == var_type
#     if np.any(idx):
#         var_error = merged["obs_error"].values[idx]
#         obs_variance[var_type - 1] = np.mean(var_error**2)
obs_variance = np.ones(19, dtype="f4")
merged["obs_variance"] = xr.DataArray(obs_variance, dims="state_variable")

# 저장
with Dataset(merged_file, "w", format="NETCDF4") as nc:
    nc.createDimension("datum", merged.dims["datum"])
    nc.createDimension("survey", merged.dims["survey"])
    nc.createDimension("state_variable", merged.dims["state_variable"])

    for var in ["obs_time", "obs_value", "obs_error", "obs_type","obs_provenance" ,
                "obs_lat", "obs_lon", "obs_depth", "obs_Xgrid", "obs_Ygrid", "obs_Zgrid"]:
        if var in merged:
            dtype = "f4" if merged[var].dtype.kind == "f" else "i4"
            v = nc.createVariable(var, dtype, ("datum",))
            v[:] = merged[var].values.astype(dtype)

    nc.createVariable("survey_time", "f4", ("survey",))[:] = merged["survey_time"].values.astype("f4")
    nc.createVariable("Nobs", "i4", ("survey",))[:] = merged["Nobs"].values.astype("i4")
    nc.createVariable("obs_variance", "f4", ("state_variable",))[:] = merged["obs_variance"].values.astype("f4")

print(f"✅ 병합 완료: {merged_file}")
