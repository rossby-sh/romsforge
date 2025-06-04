# io_utils.py

from netCDF4 import Dataset, MFDataset, num2date, date2num
from datetime import datetime

def is_netcdf4(file_path):
    with open(file_path, 'rb') as f:
        header = f.read(8)
    return header.startswith(b'CDF\x02') or header.startswith(b'\x89HDF')

def determine_open_mode(file_list):
    if isinstance(file_list, str):
        file_list = [file_list]
    return 'single' if all(is_netcdf4(f) for f in file_list) else 'mf'


def parse_time_range(date_input):
    if isinstance(date_input, str):
        t0 = datetime.fromisoformat(date_input)
        return (t0, t0)
    elif isinstance(date_input, (list, tuple)) and len(date_input) == 2:
        t_start = datetime.fromisoformat(date_input[0])
        t_end   = datetime.fromisoformat(date_input[1])
        return (t_start, t_end)
    else:
        raise ValueError("Invalid date input")

def collect_time_info(input_files, time_var, date_input, units):
    if isinstance(input_files, str):
        # initial: 단일 파일
        with Dataset(input_files) as nc:
            times = nc.variables[time_var][:]
            tdates = num2date(times, units)
            target_date = datetime.fromisoformat(date_input)
            for i, t in enumerate(tdates):
                if t == target_date:
                    return [(input_files, i, t)]
            raise ValueError("No matching date in the input file for initdate")

    elif isinstance(input_files, (list, tuple)):
        # boundary: 여러 파일
        t_start, t_end = parse_time_range(date_input)
        time_info = []

        open_mode = determine_open_mode(input_files)
        if open_mode == 'single':
            for f in input_files:
                with Dataset(f) as nc:
                    times = nc.variables[time_var][:]
                    tdates = num2date(times, units)
                    for i, t in enumerate(tdates):
                        if t_start <= t <= t_end:
                            time_info.append((f, i, t))
        elif open_mode == 'mf':
            nc = MFDataset(input_files)
            times = nc.variables[time_var][:]
            tdates = num2date(times, units)
            for i, t in enumerate(tdates):
                if t_start <= t <= t_end:
                    time_info.append((None, i, t))
        else:
            raise RuntimeError("Unknown open mode")

        return time_info

    else:
        raise ValueError("Invalid input file type: must be str or list")

