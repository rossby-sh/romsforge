# -*- coding: utf-8 -*-
"""
Created on Fri May  9 14:09:34 2025

@author: ust21
"""

import xarray as xr
import numpy as np
from netCDF4 import Dataset

pth = 'D:/shjo/ROMS_inputs/obs/pro/'
pthG = 'D:/shjo/ROMS_inputs/'

merged_file = pth + "ROMS_obs_30km.nc"


ncG = Dataset(pthG + 'roms_grd_fennel_15km_smooth_v2.nc')
lon_rho = ncG['lon_rho'][:]
lat_rho = ncG['lat_rho'][:]


NC=Dataset(merged_file)


obs_Xgrid=NC['obs_Xgrid'][:]
obs_Ygrid=NC['obs_Ygrid'][:]



np.where((obs_Xgrid < 1) | (obs_Xgrid > lon_rho.shape[1]-2))
np.any((obs_Ygrid < 1) | (obs_Ygrid > lon_rho.shape[0]-2))


lon_rho[ int(obs_Ygrid[1738977]), int(obs_Xgrid[1738977]) ]
def inspect_obs_rejection(obs_Xgrid, obs_Ygrid, obs_Zgrid, obs_type, obs_error, lon_rho, mask_rho=None):
    results = {}

    # 1. Grid 범위 벗어남
    out_x = (obs_Xgrid < 1) | (obs_Xgrid > lon_rho.shape[1] - 2)
    out_y = (obs_Ygrid < 1) | (obs_Ygrid > lon_rho.shape[0] - 2)
    results['grid_bounds_violation'] = np.sum(out_x | out_y)

    # 2. obs_error <= 0
    error_zero_or_neg = obs_error <= 0
    results['zero_or_negative_error'] = np.sum(error_zero_or_neg)

    # 3. obs_type 값 이상
    invalid_type = ~np.isin(obs_type, np.arange(1, 20))  # 1~19 사이의 값만 허용
    results['invalid_obs_type'] = np.sum(invalid_type)

    # 4. Zgrid가 0보다 작은 경우
    bad_zgrid = obs_Zgrid < 0
    results['negative_Zgrid'] = np.sum(bad_zgrid)

    # 5. 육지인지 확인 (mask_rho가 주어진 경우만)
    if mask_rho is not None:
        in_domain = (obs_Xgrid >= 0) & (obs_Xgrid < mask_rho.shape[1]) & \
                    (obs_Ygrid >= 0) & (obs_Ygrid < mask_rho.shape[0])
        mask_val = np.zeros_like(obs_Xgrid)
        mask_val[in_domain] = mask_rho[obs_Ygrid[in_domain], obs_Xgrid[in_domain]]
        on_land = mask_val == 0
        results['on_land_mask'] = np.sum(on_land)
    
    # 전체
    total = len(obs_Xgrid)
    results['total'] = total

    return results



results = inspect_obs_rejection(
    NC['obs_Xgrid'][:].astype(int),
    NC['obs_Ygrid'][:].astype(int),
    NC['obs_Zgrid'][:],
    NC['obs_type'][:],
    NC['obs_error'][:],
    lon_rho,
    mask_rho=ncG['mask_rho'][:]  # 이게 있으면 더 좋음
)

for k, v in results.items():
    print(f"{k:25}: {v}")





















