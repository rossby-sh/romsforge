import numpy as np
import datetime as dt
from netCDF4 import date2num


# Convert time to ROMS reference
def convert_time_to_ref(times, input_unit, ref_unit):
    # ref_unit = 'days since 2000-1-1 00:00:00'
    anchor = dt.datetime(2025, 1, 1, 0)
    offset = (
        date2num(anchor, ref_unit) -
        date2num(anchor, input_unit) / 86400
    )
    return times[:] / 86400 + offset

