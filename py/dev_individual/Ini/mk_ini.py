import sys
import os
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date, date2num
from py.dev_individual.libs.tools import crop_to_model_domain, load_ogcm_metadata
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import create_I as cn
import tools as tl



cfg = tl.parse_config("./config.yaml")

grd = tl.load_roms_grid(cfg.grdname)
ogcm = tl.load_ogcm_metadata(cfg.ogcm_name, cfg.ogcm_var_name)

idx, idy = tl.crop_to_model_domain(ogcm.lat, ogcm.lon, grd.lat, grd.lon)

relative_time = tl.compute_relative_time(ogcm.time[:], cfg.ref_time, ogcm.time.units)

cn.createI(cfg, grd, relative_time, ncFormat='NETCDF3_64BIT', bio_model="Fennel")
























