####################################################################################################
# Author: shjo                                                                                     #
# Description: This script merges ERA5 3-hourly forcing variables into a single NetCDF file.       #
#              Variables included: Uwind, Vwind, Tair, Qair, Cloud, sst, dqdsst, srf, lwrad,        #
#              lwrad_down, rain, Pair. Time is converted to ROMS reference time format.            #
#              Output: ERA5_3hourly_combined_2025.nc                                                #
# Date: 2025-05-26                                                                                   #
# Contact: -                                                                                         #
####################################################################################################

import sys
import os
import numpy as np
import datetime as dt
import yaml
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import create_F as cn
from utils import compute_relative_time

# Load config
with open('config_era5.yml', 'r') as f:
    config = yaml.safe_load(f)

# Construct file paths
R_FILE = os.path.join(config['input_file'])
W_FILE = os.path.join(config['output_file'])

# Load ERA5 data
ds = Dataset(R_FILE)
LON = ds.variables['longitude'][:]
LAT = np.flip(ds.variables['latitude'][:], axis=0)
TIMES = ds.variables['valid_time']

print('=== Calc ROMS Frc ===')
# Read and process variables
srf_values       = np.flip(ds.variables['ssr'][:], axis=1) / 3600
lrf_values       = np.flip(ds.variables['str'][:], axis=1) / 3600
lrf_down_values  = np.flip(ds.variables['strd'][:], axis=1) / 3600

u_values         = np.flip(ds.variables['u10'][:], axis=1)
v_values         = np.flip(ds.variables['v10'][:], axis=1)
cloud_values     = np.flip(ds.variables['tcc'][:], axis=1)
rain_values      = np.flip(ds.variables['tp'][:], axis=1) * 1000 / 3600
sst_values       = np.flip(ds.variables['sst'][:].data, axis=1) - 273.15
# sst_values = np.nan_to_num(sst_values, nan=0.0)  # NaN 제거
sst_values[sst_values != sst_values] = 0



T2_values        = np.flip(ds.variables['t2m'][:], axis=1) - 273.15
pair_values      = np.flip(ds.variables['msl'][:], axis=1) / 100

D2_values        = np.flip(ds.variables['d2m'][:], axis=1) - 273.15
TK_values        = np.flip(ds.variables['t2m'][:], axis=1)

# Calculate Qair
qair_values = 100 * (
    np.exp((17.625 * D2_values) / (243.04 + D2_values)) /
    np.exp((17.625 * T2_values) / (243.04 + T2_values))
)

# Calculate dqdsst
wind_speed = np.sqrt(u_values ** 2 + v_values ** 2)
dqdsst_values = cn.get_dqdsst(
    sst_values,
    TK_values,
    1.2,
    wind_speed,
    qair_values / 1000
)
MY_TIME_REF = config['ref_time']
# Convert time
TIME_CONVERTED_NUM = compute_relative_time(TIMES[:], TIMES.units, MY_TIME_REF)

print("=== Convert time ===")
print(num2date(TIME_CONVERTED_NUM[:5],MY_TIME_REF))
print(num2date(TIME_CONVERTED_NUM[-5:],MY_TIME_REF))

# Assemble output variables
forcing_vars = {
    'Uwind':       u_values,
    'Vwind':       v_values,
    'Tair':        T2_values,
    'Qair':        qair_values,
    'Cloud':       cloud_values,
    'sst':         sst_values,
    'dqdsst':      dqdsst_values,
    'srf':         srf_values,
    'lwrad':       lrf_values,
    'lwrad_down':  lrf_down_values,
    'rain':        rain_values,
    'Pair':        pair_values,
}

print("=== Save to netcdf ===")

# Write output NetCDF
cn.createF_era5(
    W_FILE,
    LON,
    LAT,
    TIME_CONVERTED_NUM,
    MY_TIME_REF,
    forcing_vars,
    "NETCDF3_64BIT_OFFSET"
)
#cn.createF_era5_(
#    W_FILE,
#    LON,
#    LAT,
#    TIME_CONVERTED_NUM,
#    MY_TIME_REF,
#    forcing_vars
#)
