###
#
# NFO
#
###
import sys
import os
import numpy as np
import datetime as dt
import yaml
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import ROMS_utils01 as ru  # stretching 함수가 포함된 모듈

def create_bry(cfg, grd, bry_time, bio_model=None, ncFormat='NETCDF3_CLASSIC'):
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask_rho == 1])
    if vtransform == 1 and tcline > hmin_:
        print("--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    directions = ['north', 'south', 'east', 'west']
    grids = {
        'north': 'xi', 'south': 'xi',
        'east': 'eta', 'west': 'eta'
    }

    # Dimensions
    dimensions = {
        'xi_u': L,
        'xi_v': Lp,
        'xi_rho': Lp,
        'eta_u': Mp,
        'eta_v': M,
        'eta_rho': Mp,
        's_rho': N,
        's_w': Np,
        'one': 1,
        'bry_time': len(bry_time),
        'zeta_time': len(bry_time),
        'temp_time': len(bry_time),
        'salt_time': len(bry_time),
        'v2d_time': len(bry_time),
        'v3d_time': len(bry_time)
    }

    # --- [01] Base variables + stretching 통합 정의 ---
    sc_r, Cs_r = ru.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ru.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    base_variables = {
        'spherical':    ('S1',   ('one',),   {}, 'T'),
        'Vtransform':   ('f4',   ('one',),   {'long_name': 'vertical terrain-following transformation equation'}, vtransform),
        'Vstretching':  ('f4',   ('one',),   {'long_name': 'vertical terrain-following stretching function'}, vstretching),
        'theta_s':      ('f4',   ('one',),   {'long_name': 'S-coordinate surface control parameter', 'units': 'nondimensional'}, theta_s),
        'theta_b':      ('f4',   ('one',),   {'long_name': 'S-coordinate bottom control parameter', 'units': 'nondimensional'}, theta_b),
        'Tcline':       ('f4',   ('one',),   {'long_name': 'S-coordinate surface/bottom layer width', 'units': 'meter'}, tcline),
        'hc':           ('f4',   ('one',),   {'long_name': 'S-coordinate parameter, critical depth', 'units': 'meter'}, tcline),
        'sc_r':         ('f4',   ('s_rho',), {'long_name': 'S-coordinate at RHO-points', 'units': 'nondimensional'}, sc_r),
        'Cs_r':         ('f4',   ('s_rho',), {'long_name': 'S-coordinate stretching curves at RHO-points', 'units': 'nondimensional'}, Cs_r),
        'sc_w':         ('f4',   ('s_w',),   {'long_name': 'S-coordinate at W-points', 'units': 'nondimensional'}, sc_w),
        'Cs_w':         ('f4',   ('s_w',),   {'long_name': 'S-coordinate stretching curves at W-points', 'units': 'nondimensional'}, Cs_w),
    }


    time_dims = ['bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']

    standard_vars = {
        'zeta':  ('f4', ('zeta_time', 'eta_rho', 'xi_rho'),     {'long_name': 'free-surface', 'units': 'meter'}),
        'ubar':  ('f4', ('v2d_time', 'eta_u', 'xi_u'),          {'long_name': 'vertically integrated u-momentum', 'units': 'meter second-1'}),
        'vbar':  ('f4', ('v2d_time', 'eta_v', 'xi_v'),          {'long_name': 'vertically integrated v-momentum', 'units': 'meter second-1'}),
        'u':     ('f4', ('v3d_time', 's_rho', 'eta_u', 'xi_u'), {'long_name': 'u-momentum', 'units': 'meter second-1'}),
        'v':     ('f4', ('v3d_time', 's_rho', 'eta_v', 'xi_v'), {'long_name': 'v-momentum', 'units': 'meter second-1'}),
        'temp':  ('f4', ('temp_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'potential temperature', 'units': 'Celsius'}),
        'salt':  ('f4', ('salt_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'salinity', 'units': 'PSU'})
    }

    bio_vars = {
        'npzd': {
            'NO3':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'nitrate concentration', 'units': 'mmol N m-3'}),
            'phyt':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'phytoplankton biomass', 'units': 'mmol N m-3'}),
            'zoop':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'zooplankton biomass', 'units': 'mmol N m-3'}),
            'detritus': ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'detritus concentration', 'units': 'mmol N m-3'})
        },
        'fennel': {
            'NO3':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'nitrate concentration', 'units': 'mmol N m-3'}),
            'NH4':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'ammonium concentration', 'units': 'mmol N m-3'}),
            'PO4':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'phosphate concentration', 'units': 'mmol P m-3'}),
            'chlo':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'chlorophyll concentration', 'units': 'mg Chl m-3'}),
            'phyt':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small phytoplankton biomass', 'units': 'mmol N m-3'}),
            'zoop':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'zooplankton biomass', 'units': 'mmol N m-3'}),
            'oxygen':('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'oxygen concentration', 'units': 'mmol O2 m-3'}),
            'TIC':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'total inorganic carbon', 'units': 'mmol C m-3'}),
            'alkalinity': ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'total alkalinity', 'units': 'umol kg-1'}),
            'SdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small carbon detritus', 'units': 'mmol C m-3'}),
            'LdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'large carbon detritus', 'units': 'mmol C m-3'}),
            'RdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'river carbon detritus', 'units': 'mmol C m-3'}),
            'SdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small nitrogen detritus', 'units': 'mmol N m-3'}),
            'LdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'large nitrogen detritus', 'units': 'mmol N m-3'}),
            'RdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'river nitrogen detritus', 'units': 'mmol N m-3'})
        }
    }

    # --- [02] 선택된 변수 조합 ---
    all_vars = dict(standard_vars)
    if bio_model is not None and bio_model.lower() in bio_vars:
        all_vars.update(bio_vars[bio_model.lower()])
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")
    else:
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")


    # --- [NEW] NetCDF 파일 열기 ---
    mode = 'w' if cfg.force_write else 'x'
    try:
        ncfile = Dataset(cfg.bryname, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {cfg.ininame} already exists and force_write=False ---")
        return 1

    # --- [NEW] 차원 정의 ---
    for name, size in dimensions.items():
        ncfile.createDimension(name, size)

    # --- [NEW] base_variables 기록 ---
    for name, (dtype, dims, attrs, value) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = value

    # --- [03] 방향별로 변수 확장 ---


    def get_bry_dims(dims, direction):
        keep_dim = 'xi_' if direction in ['north', 'south'] else 'eta_'
        return tuple(
            d for d in dims
            if not d.startswith(('xi_', 'eta_')) or d.startswith(keep_dim)
        )

    bry_variables = {}
    for varname, (dtype, dims, attrs) in all_vars.items():
        for d in directions:
            dims_d = get_bry_dims(dims, d)
            bry_variables[f'{varname}_{d}'] = (
                dtype, dims_d, {**attrs, 'long_name': f"{attrs['long_name']} at {d} boundary"}
            )



    # Time variables
    for name in time_dims:
        v = ncfile.createVariable(name, 'f4', (name,))
        v.units = cfg.time_ref
        v.long_name = f"time for {name.replace('_time','')} condition"
        v[:] = bry_time

    # Boundary variables
    for name, (dtype, dims, attrs) in bry_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = 0.0

    # Global attrs
    ncfile.title = cfg.global_attrs.title
    ncfile.clim_file = cfg.bryname
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs.type
    ncfile.history = f"Created on {dt.datetime.now().isoformat()}"

    ncfile.close()
    print(f"--- [+] Boundary file created: {cfg.bryname} ---")
    return 0



def create_bry_clm(cfg, grd, bry_time, bio_model=None, ncFormat='NETCDF3_CLASSIC'):
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask == 1])
    if vtransform == 1 and tcline > hmin_:
        print("--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    directions = ['north', 'south', 'east', 'west']
    grids = {
        'north': 'xi', 'south': 'xi',
        'east': 'eta', 'west': 'eta'
    }

    # Dimensions
    dimensions = {
        'xi_u': L,
        'xi_v': Lp,
        'xi_rho': Lp,
        'eta_u': Mp,
        'eta_v': M,
        'eta_rho': Mp,
        's_rho': N,
        's_w': Np,
        'one': 1,
        'bry_time': len(bry_time),
        'zeta_time': len(bry_time),
        'temp_time': len(bry_time),
        'salt_time': len(bry_time),
        'v2d_time': len(bry_time),
        'v3d_time': len(bry_time)
    }

    # --- [01] Base variables + stretching 통합 정의 ---
    sc_r, Cs_r = ru.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ru.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    base_variables = {
        'spherical':    ('S1',   ('one',),   {}, 'T'),
        'Vtransform':   ('f4',   ('one',),   {'long_name': 'vertical terrain-following transformation equation'}, vtransform),
        'Vstretching':  ('f4',   ('one',),   {'long_name': 'vertical terrain-following stretching function'}, vstretching),
        'theta_s':      ('f4',   ('one',),   {'long_name': 'S-coordinate surface control parameter', 'units': 'nondimensional'}, theta_s),
        'theta_b':      ('f4',   ('one',),   {'long_name': 'S-coordinate bottom control parameter', 'units': 'nondimensional'}, theta_b),
        'Tcline':       ('f4',   ('one',),   {'long_name': 'S-coordinate surface/bottom layer width', 'units': 'meter'}, tcline),
        'hc':           ('f4',   ('one',),   {'long_name': 'S-coordinate parameter, critical depth', 'units': 'meter'}, tcline),
        'sc_r':         ('f4',   ('s_rho',), {'long_name': 'S-coordinate at RHO-points', 'units': 'nondimensional'}, sc_r),
        'Cs_r':         ('f4',   ('s_rho',), {'long_name': 'S-coordinate stretching curves at RHO-points', 'units': 'nondimensional'}, Cs_r),
        'sc_w':         ('f4',   ('s_w',),   {'long_name': 'S-coordinate at W-points', 'units': 'nondimensional'}, sc_w),
        'Cs_w':         ('f4',   ('s_w',),   {'long_name': 'S-coordinate stretching curves at W-points', 'units': 'nondimensional'}, Cs_w),
    }


    time_dims = ['bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']

    standard_vars = {
        'zeta':  ('f4', ('zeta_time', 'eta_rho', 'xi_rho'),     {'long_name': 'free-surface', 'units': 'meter'}),
        'ubar':  ('f4', ('v2d_time', 'eta_u', 'xi_u'),          {'long_name': 'vertically integrated u-momentum', 'units': 'meter second-1'}),
        'vbar':  ('f4', ('v2d_time', 'eta_v', 'xi_v'),          {'long_name': 'vertically integrated v-momentum', 'units': 'meter second-1'}),
        'u':     ('f4', ('v3d_time', 's_rho', 'eta_u', 'xi_u'), {'long_name': 'u-momentum', 'units': 'meter second-1'}),
        'v':     ('f4', ('v3d_time', 's_rho', 'eta_v', 'xi_v'), {'long_name': 'v-momentum', 'units': 'meter second-1'}),
        'temp':  ('f4', ('temp_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'potential temperature', 'units': 'Celsius'}),
        'salt':  ('f4', ('salt_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'salinity', 'units': 'PSU'})
    }

    bio_vars = {
        'npzd': {
            'NO3':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'nitrate concentration', 'units': 'mmol N m-3'}),
            'phyt':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'phytoplankton biomass', 'units': 'mmol N m-3'}),
            'zoop':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'zooplankton biomass', 'units': 'mmol N m-3'}),
            'detritus': ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'detritus concentration', 'units': 'mmol N m-3'})
        },
        'fennel': {
            'NO3':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'nitrate concentration', 'units': 'mmol N m-3'}),
            'NH4':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'ammonium concentration', 'units': 'mmol N m-3'}),
            'PO4':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'phosphate concentration', 'units': 'mmol P m-3'}),
            'chlo':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'chlorophyll concentration', 'units': 'mg Chl m-3'}),
            'phyt':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small phytoplankton biomass', 'units': 'mmol N m-3'}),
            'zoop':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'zooplankton biomass', 'units': 'mmol N m-3'}),
            'oxygen':('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'oxygen concentration', 'units': 'mmol O2 m-3'}),
            'TIC':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'total inorganic carbon', 'units': 'mmol C m-3'}),
            'alkalinity': ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'total alkalinity', 'units': 'umol kg-1'}),
            'SdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small carbon detritus', 'units': 'mmol C m-3'}),
            'LdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'large carbon detritus', 'units': 'mmol C m-3'}),
            'RdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'river carbon detritus', 'units': 'mmol C m-3'}),
            'SdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small nitrogen detritus', 'units': 'mmol N m-3'}),
            'LdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'large nitrogen detritus', 'units': 'mmol N m-3'}),
            'RdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'river nitrogen detritus', 'units': 'mmol N m-3'})
        }
    }

    # --- [02] 선택된 변수 조합 ---
    all_vars = dict(standard_vars)
    if bio_model is not None and bio_model.lower() in bio_vars:
        all_vars.update(bio_vars[bio_model.lower()])
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")
    else:
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")


    # --- [NEW] NetCDF 파일 열기 ---
    mode = 'w' if cfg.force_write else 'x'
    try:
        ncfile = Dataset(cfg.bryname, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {cfg.ininame} already exists and force_write=False ---")
        return 1

    # --- [NEW] 차원 정의 ---
    for name, size in dimensions.items():
        ncfile.createDimension(name, size)

    # --- [NEW] base_variables 기록 ---
    for name, (dtype, dims, attrs, value) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = value

    # --- [03] 방향별로 변수 확장 ---

    def get_bry_dims(dims, direction):
        keep_dim = 'xi_' if direction in ['north', 'south'] else 'eta_'
        return tuple(
            d for d in dims
            if not d.startswith(('xi_', 'eta_')) or d.startswith(keep_dim)
        )

    bry_variables = {}
    for varname, (dtype, dims, attrs) in all_vars.items():
        for d in directions:
            dims_d = get_bry_dims(dims, d)
            bry_variables[f'{varname}_{d}'] = (
                dtype, dims_d, {**attrs, 'long_name': f"{attrs['long_name']} at {d} boundary"}
            )



    # Time variables
    for name in time_dims:
        v = ncfile.createVariable(name, 'f4', (name,))
        v.units = "day"
        v.cycle_length = 365.25
        v.long_name = f"time for {name.replace('_time','')} condition"
        v[:] = bry_time

    # Boundary variables
    for name, (dtype, dims, attrs) in bry_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = 0.0

    # Global attrs
    ncfile.title = cfg.global_attrs.title
    ncfile.clim_file = cfg.bryname
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs.type
    ncfile.history = f"Created on {dt.datetime.now().isoformat()}"

    ncfile.close()
    print(f"--- [+] Boundary file created: {cfg.bryname} ---")
    return 0


def create_bry_nifs(bryname,cfg, grd, bry_time, bio_model=None, ncFormat='NETCDF3_CLASSIC'):
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask_rho == 1])
    if vtransform == 1 and tcline > hmin_:
        print("--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    directions = ['north', 'south', 'east', 'west']
    grids = {
        'north': 'xi', 'south': 'xi',
        'east': 'eta', 'west': 'eta'
    }

    # Dimensions
    dimensions = {
        'xi_u': L,
        'xi_v': Lp,
        'xi_rho': Lp,
        'eta_u': Mp,
        'eta_v': M,
        'eta_rho': Mp,
        's_rho': N,
        's_w': Np,
        'one': 1,
        'bry_time': len(bry_time),
        'zeta_time': len(bry_time),
        'temp_time': len(bry_time),
        'salt_time': len(bry_time),
        'v2d_time': len(bry_time),
        'v3d_time': len(bry_time)
    }

    # --- [01] Base variables + stretching 통합 정의 ---
    sc_r, Cs_r = ru.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ru.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    base_variables = {
        'spherical':    ('S1',   ('one',),   {}, 'T'),
        'Vtransform':   ('f4',   ('one',),   {'long_name': 'vertical terrain-following transformation equation'}, vtransform),
        'Vstretching':  ('f4',   ('one',),   {'long_name': 'vertical terrain-following stretching function'}, vstretching),
        'theta_s':      ('f4',   ('one',),   {'long_name': 'S-coordinate surface control parameter', 'units': 'nondimensional'}, theta_s),
        'theta_b':      ('f4',   ('one',),   {'long_name': 'S-coordinate bottom control parameter', 'units': 'nondimensional'}, theta_b),
        'Tcline':       ('f4',   ('one',),   {'long_name': 'S-coordinate surface/bottom layer width', 'units': 'meter'}, tcline),
        'hc':           ('f4',   ('one',),   {'long_name': 'S-coordinate parameter, critical depth', 'units': 'meter'}, tcline),
        'sc_r':         ('f4',   ('s_rho',), {'long_name': 'S-coordinate at RHO-points', 'units': 'nondimensional'}, sc_r),
        'Cs_r':         ('f4',   ('s_rho',), {'long_name': 'S-coordinate stretching curves at RHO-points', 'units': 'nondimensional'}, Cs_r),
        'sc_w':         ('f4',   ('s_w',),   {'long_name': 'S-coordinate at W-points', 'units': 'nondimensional'}, sc_w),
        'Cs_w':         ('f4',   ('s_w',),   {'long_name': 'S-coordinate stretching curves at W-points', 'units': 'nondimensional'}, Cs_w),
    }


    time_dims = ['bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']

    standard_vars = {
        'zeta':  ('f4', ('zeta_time', 'eta_rho', 'xi_rho'),     {'long_name': 'free-surface', 'units': 'meter'}),
        'ubar':  ('f4', ('v2d_time', 'eta_u', 'xi_u'),          {'long_name': 'vertically integrated u-momentum', 'units': 'meter second-1'}),
        'vbar':  ('f4', ('v2d_time', 'eta_v', 'xi_v'),          {'long_name': 'vertically integrated v-momentum', 'units': 'meter second-1'}),
        'u':     ('f4', ('v3d_time', 's_rho', 'eta_u', 'xi_u'), {'long_name': 'u-momentum', 'units': 'meter second-1'}),
        'v':     ('f4', ('v3d_time', 's_rho', 'eta_v', 'xi_v'), {'long_name': 'v-momentum', 'units': 'meter second-1'}),
        'temp':  ('f4', ('temp_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'potential temperature', 'units': 'Celsius'}),
        'salt':  ('f4', ('salt_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'salinity', 'units': 'PSU'})
    }

    bio_vars = {
        'npzd': {
            'NO3':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'nitrate concentration', 'units': 'mmol N m-3'}),
            'phyt':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'phytoplankton biomass', 'units': 'mmol N m-3'}),
            'zoop':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'zooplankton biomass', 'units': 'mmol N m-3'}),
            'detritus': ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'detritus concentration', 'units': 'mmol N m-3'})
        },
        'fennel': {
            'NO3':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'nitrate concentration', 'units': 'mmol N m-3'}),
            'NH4':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'ammonium concentration', 'units': 'mmol N m-3'}),
            'PO4':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'phosphate concentration', 'units': 'mmol P m-3'}),
            'chlo':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'chlorophyll concentration', 'units': 'mg Chl m-3'}),
            'phyt':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small phytoplankton biomass', 'units': 'mmol N m-3'}),
            'zoop':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'zooplankton biomass', 'units': 'mmol N m-3'}),
            'oxygen':('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'oxygen concentration', 'units': 'mmol O2 m-3'}),
            'TIC':   ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'total inorganic carbon', 'units': 'mmol C m-3'}),
            'alkalinity': ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'total alkalinity', 'units': 'umol kg-1'}),
            'SdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small carbon detritus', 'units': 'mmol C m-3'}),
            'LdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'large carbon detritus', 'units': 'mmol C m-3'}),
            'RdeC':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'river carbon detritus', 'units': 'mmol C m-3'}),
            'SdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'small nitrogen detritus', 'units': 'mmol N m-3'}),
            'LdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'large nitrogen detritus', 'units': 'mmol N m-3'}),
            'RdeN':  ('f4', ('bry_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'river nitrogen detritus', 'units': 'mmol N m-3'})
        }
    }

    # --- [02] 선택된 변수 조합 ---
    all_vars = dict(standard_vars)
    if bio_model is not None and bio_model.lower() in bio_vars:
        all_vars.update(bio_vars[bio_model.lower()])
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")
    else:
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")


    # --- [NEW] NetCDF 파일 열기 ---
    mode = 'w' if cfg.force_write else 'x'
    try:
        ncfile = Dataset(bryname, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {bryname} already exists and force_write=False ---")
        return 1

    # --- [NEW] 차원 정의 ---
    for name, size in dimensions.items():
        ncfile.createDimension(name, size)

    # --- [NEW] base_variables 기록 ---
    for name, (dtype, dims, attrs, value) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = value

    # --- [03] 방향별로 변수 확장 ---


    def get_bry_dims(dims, direction):
        keep_dim = 'xi_' if direction in ['north', 'south'] else 'eta_'
        return tuple(
            d for d in dims
            if not d.startswith(('xi_', 'eta_')) or d.startswith(keep_dim)
        )

    bry_variables = {}
    for varname, (dtype, dims, attrs) in all_vars.items():
        for d in directions:
            dims_d = get_bry_dims(dims, d)
            bry_variables[f'{varname}_{d}'] = (
                dtype, dims_d, {**attrs, 'long_name': f"{attrs['long_name']} at {d} boundary"}
            )



    # Time variables
    for name in time_dims:
        v = ncfile.createVariable(name, 'f4', (name,))
        v.units = cfg.time_ref
        v.long_name = f"time for {name.replace('_time','')} condition"
        v[:] = bry_time

    # Boundary variables
    for name, (dtype, dims, attrs) in bry_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = 0.0

    # Global attrs
    ncfile.title = cfg.global_attrs_bry.title
    ncfile.clim_file = bryname
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs_bry.type
    ncfile.history = f"Created on {dt.datetime.now().isoformat()}"

    ncfile.close()
    print(f"--- [+] Boundary file created: {bryname} ---")
    return 0







def create_bry_tmp(cfg, grd, bry_time, bio_model=None, ncFormat='NETCDF3_CLASSIC'):
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask == 1])
    if vtransform == 1 and tcline > hmin_:
        print("--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    directions = ['north', 'south', 'east', 'west']
    grids = {
        'north': 'xi', 'south': 'xi',
        'east': 'eta', 'west': 'eta'
    }

    # Dimensions
    dimensions = {
        'xi_u': L,
        'xi_v': Lp,
        'xi_rho': Lp,
        'eta_u': Mp,
        'eta_v': M,
        'eta_rho': Mp,
        's_rho': N,
        's_w': Np,
        'one': 1,
        'bry_time': len(bry_time),
        'zeta_time': len(bry_time),
        'temp_time': len(bry_time),
        'salt_time': len(bry_time),
        'v2d_time': len(bry_time),
        'v3d_time': len(bry_time)
    }

    # Base variables (1D control variables)
    base_variables = {
        'spherical':    ('S1',   ('one',),   {}, 'T'),
        'Vtransform':   ('f4',   ('one',),   {'long_name': 'vertical terrain-following transformation equation'}, vtransform),
        'Vstretching':  ('f4',   ('one',),   {'long_name': 'vertical terrain-following stretching function'}, vstretching),
        'theta_s':      ('f4',   ('one',),   {'long_name': 'S-coordinate surface control parameter', 'units': 'nondimensional'}, theta_s),
        'theta_b':      ('f4',   ('one',),   {'long_name': 'S-coordinate bottom control parameter', 'units': 'nondimensional'}, theta_b),
        'Tcline':       ('f4',   ('one',),   {'long_name': 'S-coordinate surface/bottom layer width', 'units': 'meter'}, tcline),
        'hc':           ('f4',   ('one',),   {'long_name': 'S-coordinate parameter, critical depth', 'units': 'meter'}, tcline),
    }

    # Stretching
    sc_r, Cs_r = ru.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ru.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    # Stretching variables
    base_variables.update({
        'sc_r': ('f4', ('s_rho',), {'long_name': 'S-coordinate at RHO-points', 'units': 'nondimensional'}, sc_r),
        'Cs_r': ('f4', ('s_rho',), {'long_name': 'S-coordinate stretching curves at RHO-points', 'units': 'nondimensional'}, Cs_r),
        'sc_w': ('f4', ('s_w',),   {'long_name': 'S-coordinate at W-points', 'units': 'nondimensional'}, sc_w),
        'Cs_w': ('f4', ('s_w',),   {'long_name': 'S-coordinate stretching curves at W-points', 'units': 'nondimensional'}, Cs_w),
    })

    time_dims = ['bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']

    bry_variables = {}
    for d in directions:
        g = grids[d]
        bry_variables.update({
            f'zeta_{d}': ('zeta_time', f'{g}_rho'),
            f'ubar_{d}': ('v2d_time', f'{g}_u'),
            f'vbar_{d}': ('v2d_time', f'{g}_v'),
            f'temp_{d}': ('temp_time', 's_rho', f'{g}_rho'),
            f'salt_{d}': ('salt_time', 's_rho', f'{g}_rho'),
            f'u_{d}':    ('v3d_time', 's_rho', f'{g}_u'),
            f'v_{d}':    ('v3d_time', 's_rho', f'{g}_v')
        })

    # Biological variables
    bio_tracers = {
        'npzd': {
            'NO3': 'nitrate concentration',
            'phyt': 'phytoplankton biomass',
            'zoop': 'zooplankton biomass',
            'detritus': 'detritus concentration'
        },
        'fennel': {
            'NO3': 'nitrate concentration',
            'NH4': 'ammonium concentration',
            'PO4': 'phosphate concentration',
            'chlo': 'chlorophyll concentration',
            'phyt': 'small phytoplankton biomass',
            'zoop': 'zooplankton biomass',
            'oxygen': 'oxygen concentration',
            'TIC': 'total inorganic carbon',
            'alkalinity': 'total alkalinity',
            'SdeC': 'small carbon-detritus concentration',
            'LdeC': 'large carbon-detritus concentration',
            'RdeC': 'river carbon-detritus concentration',
            'SdeN': 'small nitrogen-detritus concentration',
            'LdeN': 'large nitrogen-detritus concentration',
            'RdeN': 'river nitrogen-detritus concentration'
        }
    }

    if bio_model in bio_tracers:
        for tracer, desc in bio_tracers[bio_model].items():
            for d in directions:
                g = grids[d]
                bry_variables[f'{tracer}_{d}'] = ('bry_time', 's_rho', f'{g}_rho')

    # Create NetCDF
    mode = 'w' if cfg.force_write else 'x'
    try:
        ncfile = Dataset(cfg.bryname, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {cfg.bryname} already exists and force_write=False ---")
        return 1

    for name, size in dimensions.items():
        ncfile.createDimension(name, size)

    for name, (dtype, dims, attrs, value) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)
        var[:] = value

    # Time variables
    for name in time_dims:
        v = ncfile.createVariable(name, 'f4', (name,))
        v.units = cfg.time_ref
        v.long_name = f"time for {name.replace('_time','')} condition"
        v[:] = bry_time

    # Boundary variables
    for name, dims in bry_variables.items():
        var = ncfile.createVariable(name, 'f4', dims)
        var.long_name = name.replace('_', ' ')
        var.units = 'unknown'
        var[:] = 0.0

    # Global attrs
    ncfile.title = cfg.global_attrs.title
    ncfile.clim_file = cfg.bryname
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs.type
    ncfile.history = cfg.global_attrs.history

    ncfile.close()
    print(f"✅ Boundary file created: {cfg.bryname}")
    return 0


def createB(
    filename,
    topo,
    mask,
    MyVar,
    time_array,
    time_ref,
    NSEW,
    bio_model='Fennel',
    ncFormat='NETCDF3_CLASSIC'
):
    hmin_ = np.min(topo[mask == 1])
    if MyVar['Vtransform'] == 1 and MyVar['Tcline'] > hmin_:
        raise ValueError("Tcline is greater than minimum depth under Vtransform=1")

    Mp, Lp = topo.shape
    L, M, Np = Lp - 1, Mp - 1, MyVar['Layer_N'] + 1

    sc_r, Cs_r = ru.stretching(MyVar['Vstretching'], MyVar['Theta_s'], MyVar['Theta_b'], MyVar['Layer_N'], 0)
    sc_w, Cs_w = ru.stretching(MyVar['Vstretching'], MyVar['Theta_s'], MyVar['Theta_b'], MyVar['Layer_N'], 1)

    MyVar['sc_r'] = sc_r
    MyVar['sc_w'] = sc_w
    MyVar['Cs_r'] = Cs_r
    MyVar['Cs_w'] = Cs_w

    ncfile = Dataset(filename, mode='w', format=ncFormat)

    ncfile.createDimension('xi_u', L)
    ncfile.createDimension('xi_v', Lp)
    ncfile.createDimension('xi_rho', Lp)
    ncfile.createDimension('eta_u', Mp)
    ncfile.createDimension('eta_v', M)
    ncfile.createDimension('eta_rho', Mp)
    ncfile.createDimension('s_rho', MyVar['Layer_N'])
    ncfile.createDimension('s_w', Np)
    ncfile.createDimension('one', 1)
    ncfile.createDimension('bry_time', len(time_array))
    ncfile.createDimension('zeta_time', len(time_array))
    ncfile.createDimension('temp_time', len(time_array))
    ncfile.createDimension('salt_time', len(time_array))
    ncfile.createDimension('v2d_time', len(time_array))
    ncfile.createDimension('v3d_time', len(time_array))

    time_var = ncfile.createVariable('bry_time', 'f4', ('bry_time',))
    time_var.long_name = 'time for boundary condition'
    time_var.units = time_ref
    time_var[:] = time_array

    for time_name in ['zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time']:
        time_var = ncfile.createVariable(time_name, 'f4', (time_name,))
        time_var.long_name = 'time for ' + time_name.replace('_time', '') + ' condition'
        time_var.units = time_ref
        time_var[:] = time_array

    one_vars = {
        'spherical':    ('S1', 'grid type logical switch', '', 'T,F', 'spherical cartesian', 'T'),
        'Vtransform':   ('f4', 'vertical terrain-following transformation equation', '', '', '', MyVar['Vtransform']),
        'Vstretching':  ('f4', 'vertical terrain-following stretching function', '', '', '', MyVar['Vstretching']),
#        'tstart':       ('f4', 'start processing day', 'day', '', '', MyVar['tstart']),
#        'tend':         ('f4', 'end processing day', 'day', '', '', MyVar['tend']),
        'theta_s':      ('f4', 'S-coordinate surface control parameter', 'nondimensional', '', '', MyVar['Theta_s']),
        'theta_b':      ('f4', 'S-coordinate bottom control parameter', 'nondimensional', '', '', MyVar['Theta_b']),
        'Tcline':       ('f4', 'S-coordinate surface/bottom layer width', 'meter', '', '', MyVar['Tcline']),
        'hc':           ('f4', 'S-coordinate parameter, critical depth', 'meter', '', '', MyVar['Tcline'])
    }

    for varname, meta in one_vars.items():
        ncfile.createVariable(varname, meta[0], ('one',))
        ncfile[varname].long_name = meta[1]
        if meta[2]:
            ncfile[varname].units = meta[2]
        if varname == 'spherical':
            ncfile[varname].flag_values = meta[3]
            ncfile[varname].flag_meanings = meta[4]
        ncfile[varname][:] = meta[5]

    for name, arr in zip(['sc_r', 'sc_w', 'Cs_r', 'Cs_w'], [sc_r, sc_w, Cs_r, Cs_w]):
        dim = 's_rho' if 'r' in name else 's_w'
        var = ncfile.createVariable(name, 'f4', (dim,))
        var.long_name = f'S-coordinate {"stretching curves" if "Cs" in name else "at " + ("RHO" if dim == "s_rho" else "W") + "-points"}'
        if 'Cs' in name:
            var.units = 'nondimensional'
            var.valid_min = -1.
            var.valid_max = 0.
        else:
            var.valid_min = -1.
            var.valid_max = 0.
            var.positive = 'up'
            var.standard_name = 'ocena_s_coordinate_g1' if MyVar['Vtransform'] == 1 else 'ocena_s_coordinate_g2'
            var.formula_terms = f's: {dim} C: Cs_{dim[-1]} eta: zeta depth: h depth_c: hc'
        var[:] = arr

    directions = ['north', 'south', 'east', 'west']
    grids = ['xi', 'xi', 'eta', 'eta']
    for d, g, is_active in zip(directions, grids, NSEW):
        if not is_active:
            continue
        var_list = {
            f'zeta_{d}':    ('zeta_time', f'{g}_rho', None, 'free-surface', 'meter', 'lon_rho'),
            f'ubar_{d}':    ('v2d_time', f'{g}_u', None, 'vertically integrated u-momentum component', 'meter second-1', 'lon_u'),
            f'vbar_{d}':    ('v2d_time', f'{g}_v', None, 'vertically integrated v-momentum component', 'meter second-1', 'lon_v'),
            f'temp_{d}':    ('temp_time', 's_rho', f'{g}_rho', 'potential temperature', 'Celsius', 'lon_rho'),
            f'salt_{d}':    ('salt_time', 's_rho', f'{g}_rho', 'salinity', 'PSU', 'lon_rho'),
            f'u_{d}':       ('v3d_time', 's_rho', f'{g}_u', 'u-momentum component', 'meter second-1', 'lon_u'),
            f'v_{d}':       ('v3d_time', 's_rho', f'{g}_v', 'v-momentum component', 'meter second-1', 'lon_v'),
            f'NO3_{d}':     ('bry_time', 's_rho', f'{g}_rho', 'nitrate concentration', 'millimole nitrogen meter-3', 'lon_rho'),
            f'NH4_{d}':     ('bry_time', 's_rho', f'{g}_rho', 'ammonium concentration', 'millimole nitrogen meter-3', 'lon_rho'),
            f'PO4_{d}':     ('bry_time', 's_rho', f'{g}_rho', 'ammonium concentration', 'millimole po4 meter-3', 'lon_rho'),
            f'chlo_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'chlorophyll concentration', 'millimole chlorphyll meter-3', 'lon_rho'),
            f'phyt_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'small phytoplankton biomass', 'millimole nitrogen meter-3', 'lon_rho'),
            f'zoop_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'zooplankton biomass', 'millimole nitrogen meter-3', 'lon_rho'),
            f'oxygen_{d}':  ('bry_time', 's_rho', f'{g}_rho', 'oxygen concentration', 'millimole oxygen meter-3', 'lon_rho'),
            f'TIC_{d}':     ('bry_time', 's_rho', f'{g}_rho', 'total inorganic carbon', 'millimole carbon meter-3', 'lon_rho'),
            f'alkalinity_{d}': ('bry_time', 's_rho', f'{g}_rho', 'total alkalinity', 'milliequivalent meter-3', 'lon_rho'),
            f'SdeC_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'small carbon-detritus concentration', 'millimole carbon meter-3', 'lon_rho'),
            f'LdeC_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'large carbon-detritus concentration', 'millimole carbon meter-3', 'lon_rho'),
            f'RdeC_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'river carbon-detritus concentration', 'millimole carbon meter-3', 'lon_rho'),
            f'SdeN_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'small nitrogen-detritus concentration', 'millimole nitrogen meter-3', 'lon_rho'),
            f'LdeN_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'large nitrogen-detritus concentration', 'millimole nitrogen meter-3', 'lon_rho'),
            f'RdeN_{d}':    ('bry_time', 's_rho', f'{g}_rho', 'river nitrogen-detritus concentration', 'millimole nitrogen meter-3', 'lon_rho')
        }

        for vname, dims in var_list.items():
            if dims[0] not in ncfile.dimensions:
                ncfile.createDimension(dims[0], len(time_array))
            valid_dims = tuple(d for d in dims[:3] if d)
            var = ncfile.createVariable(vname, 'f4', valid_dims)
            var.long_name = dims[3]
            var.units = dims[4]
            if 's_rho' in valid_dims:
                var.coordinates = f"{dims[5]} s_rho {dims[0]}"
            else:
                var.coordinates = f"{dims[5]} {dims[0]}"
            var[:] = 0.0

    ncfile.close()
    print(f"✅ Boundary file created: {filename}")




def createB_NPZD(My_Bry, mask, topo, MyVar, NSEW, Bry_time, My_time_ref, Title, ncFormat='NETCDF3_CLASSIC'):

    hmin_ = np.min(topo[mask == 1])
    if MyVar['Vtransform'] == 1 and MyVar['Tcline'] > hmin_:
        raise

    Mp, Lp = topo.shape
    L, M, Np = Lp - 1, Mp - 1, MyVar['Layer_N'] + 1

    ncfile = Dataset(My_Bry, mode='w', format=ncFormat)

    # Dimensions
    ncfile.createDimension('xi_u', L)
    ncfile.createDimension('xi_v', Lp)
    ncfile.createDimension('xi_rho', Lp)
    ncfile.createDimension('eta_u', Mp)
    ncfile.createDimension('eta_v', M)
    ncfile.createDimension('eta_rho', Mp)
    ncfile.createDimension('s_rho', MyVar['Layer_N'])
    ncfile.createDimension('s_w', Np)
    ncfile.createDimension('tracer', 2)
    ncfile.createDimension('one', 1)

    for time_name in ["bry_time", "zeta_time", "temp_time", "salt_time", "v2d_time", "v3d_time", "NO3_time", "phyt_time", "zoop_time", "detritus_time"]:
        ncfile.createDimension(time_name, len(Bry_time))
        tvar = ncfile.createVariable(time_name, 'f4', (time_name,))
        tvar.long_name = 'Time for boundary'
        tvar.units = My_time_ref

    # Scalars
    ncfile.createVariable('spherical', 'S1', ('one',)).setncattr('long_name', 'grid type logical switch')
    ncfile['spherical'].flag_values = 'T,F'
    ncfile['spherical'].flag_meanings = 'spherical cartesian'

    for varname in ['Vtransform', 'Vstretching', 'tstart', 'tend', 'theta_s', 'theta_b', 'Tcline', 'hc']:
        v = ncfile.createVariable(varname, 'f4', ('one',))
        v.long_name = varname.replace('_', ' ')
        v.units = 'meter' if varname in ['Tcline', 'hc'] else 'nondimensional' if 'theta' in varname else 'day'

    ncfile.createVariable('sc_r', 'f4', ('s_rho',)).setncatts({
        'long_name': 'S-coordinate at RHO-points',
        'valid_min': -1., 'valid_max': 0., 'positive': 'up',
        'formula_terms': 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc',
        'standard_name': f'ocena_s_coordinate_g{MyVar["Vtransform"]}'
    })
    ncfile.createVariable('sc_w', 'f4', ('s_w',)).setncatts({
        'long_name': 'S-coordinate at W-points',
        'valid_min': -1., 'valid_max': 0., 'positive': 'up',
        'formula_terms': 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc',
        'standard_name': f'ocena_s_coordinate_g{MyVar["Vtransform"]}'
    })
    ncfile.createVariable('Cs_r', 'f4', ('s_rho',)).setncatts({
        'long_name': 'S-coordinate stretching curves at RHO-points',
        'units': 'nondimensional', 'valid_min': -1., 'valid_max': 0.
    })
    ncfile.createVariable('Cs_w', 'f4', ('s_w',)).setncatts({
        'long_name': 'S-coordinate stretching curves at W-points',
        'units': 'nondimensional', 'valid_min': -1., 'valid_max': 0.
    })

    # NSEW boundary variables
    directions = ['north', 'south', 'east', 'west']
    grids = ['xi', 'xi', 'eta', 'eta']

    for m, use in enumerate(NSEW):
        if not use:
            continue
        d = directions[m]
        g = grids[m]

        print(f'!!! Make bry: {d} !!!')

        var_list = {
            f'u_{d}': ('v3d_time', 's_rho', f'{g}_u', 'u-momentum component', 'meter second-1', 'lon_u'),
            f'v_{d}': ('v3d_time', 's_rho', f'{g}_v', 'v-momentum component', 'meter second-1', 'lon_v'),
            f'temp_{d}': ('temp_time', 's_rho', f'{g}_rho', 'potential temperature', 'Celsius', 'lon_rho'),
            f'salt_{d}': ('salt_time', 's_rho', f'{g}_rho', 'salinity', 'PSU', 'lon_rho'),
            f'ubar_{d}': ('v2d_time', f'{g}_u', '', 'vertically integrated u-momentum component', 'meter second-1', 'lon_u'),
            f'vbar_{d}': ('v2d_time', f'{g}_v', '', 'vertically integrated v-momentum component', 'meter second-1', 'lon_v'),
            f'zeta_{d}': ('zeta_time', f'{g}_rho', '', 'free-surface', 'meter', 'lon_rho'),
            f'phyt_{d}': ('phyt_time', 's_rho', f'{g}_rho', 'phytoplankton', 'millimole phytoplankton meter-3', 'lon_rho'),
            f'NO3_{d}': ('NO3_time', 's_rho', f'{g}_rho', 'NO3', 'millimole nitrogen meter-3', 'lon_rho'),
            f'zoop_{d}': ('zoop_time', 's_rho', f'{g}_rho', 'zooplankton', 'millimole_zooplankton meter-3', 'lon_rho'),
            f'detritus_{d}': ('detritus_time', 's_rho', f'{g}_rho', 'detritus', 'millimole nitrogen meter-3', 'lon_rho')
        }

        for vname, dims in var_list.items():
            dim = tuple(d for d in dims[:3] if d)
            var = ncfile.createVariable(vname, 'f4', dim)
            var.long_name = dims[3] if 'component' in dims[3] or dims[3] == 'free-surface' else f'{dims[3]} {d} boundary condition'
            var.units = dims[4]
            var.coordinates = f'{dims[5]} s_rho bry_time' if 's_rho' in dim else f'{dims[5]} bry_time'
            var[:] = 0

    # Set scalar values
    sc_r, Cs_r = ru.stretching(MyVar['Vstretching'], MyVar['Theta_s'], MyVar['Theta_b'], MyVar['Layer_N'], 0)
    sc_w, Cs_w = ru.stretching(MyVar['Vstretching'], MyVar['Theta_s'], MyVar['Theta_b'], MyVar['Layer_N'], 1)

    ncfile['spherical'][:] = 'T'
    ncfile['Vtransform'][:] = MyVar['Vtransform']
    ncfile['Vstretching'][:] = MyVar['Vstretching']
    ncfile['tstart'][:] = 0
    ncfile['tend'][:] = 0
    ncfile['theta_s'][:] = MyVar['Theta_s']
    ncfile['theta_b'][:] = MyVar['Theta_b']
    ncfile['Tcline'][:] = MyVar['Tcline']
    ncfile['hc'][:] = MyVar['Tcline']
    ncfile['sc_r'][:] = sc_r
    ncfile['Cs_r'][:] = Cs_r
    ncfile['sc_w'][:] = sc_w
    ncfile['Cs_w'][:] = Cs_w

    for t in ["bry_time", "zeta_time", "temp_time", "salt_time", "v2d_time", "v3d_time", "NO3_time", "phyt_time", "zoop_time", "detritus_time"]:
        ncfile[t][:] = Bry_time

    ncfile.title = Title
    ncfile.clim_file = My_Bry
    ncfile.grd_file = ''
    ncfile.type = 'bry file'
    ncfile.history = 'ROMS'
    ncfile.close()

