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

def create_ini(cfg, grd, initime_num, ncFormat='NETCDF3_CLASSIC', bio_model=None):


    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask == 1])
    if vtransform == 1 and tcline > hmin_:
        print(f"[ERROR](Tcline must be <= hmin when Vtransform == 1")
        return 1 


    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    dimensions = {
        'xi_u':       L,
        'xi_v':       Lp,
        'xi_rho':     Lp,
        'eta_u':      Mp,
        'eta_v':      M,
        'eta_rho':    Mp,
        's_rho':      N,
        's_w':        Np,
        'tracer':     2,
        'one':        1,
        'ocean_time': None,
    }

    base_variables = {
        'spherical':     ('S1',   ('one',),   {},),
        'Vtransform':    ('f4',   ('one',),   {'long_name': 'vertical terrain-following transformation equation'}),
        'Vstretching':   ('f4',   ('one',),   {'long_name': 'vertical terrain-following stretching function'}),
        'theta_s':       ('f4',   ('one',),   {'long_name': 'S-coordinate surface control parameter', 'units': 'nondimensional'}),
        'theta_b':       ('f4',   ('one',),   {'long_name': 'S-coordinate bottom control parameter', 'units': 'nondimensional'}),
        'Tcline':        ('f4',   ('one',),   {'long_name': 'S-coordinate surface/bottom layer width', 'units': 'meter'}),
        'hc':            ('f4',   ('one',),   {'long_name': 'S-coordinate parameter, critical depth', 'units': 'meter'}),

        'sc_r':          ('f4',   ('s_rho',), {'long_name': 'S-coordinate at RHO-points', 'units': 'nondimensional'}),
        'Cs_r':          ('f4',   ('s_rho',), {'long_name': 'S-coordinate stretching curves at RHO-points', 'units': 'nondimensional'}),
        'sc_w':          ('f4',   ('s_w',),   {'long_name': 'S-coordinate at W-points', 'units': 'nondimensional'}),
        'Cs_w':          ('f4',   ('s_w',),   {'long_name': 'S-coordinate stretching curves at W-points', 'units': 'nondimensional'}),

        'ocean_time':    ('f4',   ('ocean_time',), {'units': cfg.time_ref}),

        'u':             ('f4',   ('ocean_time', 's_rho', 'eta_u', 'xi_u'), {'long_name': 'u-momentum component', 'units': 'meter second-1'}),
        'v':             ('f4',   ('ocean_time', 's_rho', 'eta_v', 'xi_v'), {'long_name': 'v-momentum component', 'units': 'meter second-1'}),
        'ubar':          ('f4',   ('ocean_time', 'eta_u', 'xi_u'),         {'long_name': 'vertically integrated u-momentum component', 'units': 'meter second-1'}),
        'vbar':          ('f4',   ('ocean_time', 'eta_v', 'xi_v'),         {'long_name': 'vertically integrated v-momentum component', 'units': 'meter second-1'}),
        'zeta':          ('f4',   ('ocean_time', 'eta_rho', 'xi_rho'),     {'long_name': 'free-surface', 'units': 'meter'}),
        'temp':          ('f4',   ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'potential temperature', 'units': 'Celsius'}),
        'salt':          ('f4',   ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), {'long_name': 'salinity', 'units': 'PSU'}),
    }

    # 추가 tracer 변수 정의
    bio_tracers = {
        'NPZD': {
            'NO3':           {'long_name': 'nitrate concentration', 'units': 'millimole nitrogen meter-3'},
            'phytoplankton': {'long_name': 'phytoplankton biomass',  'units': 'millimole nitrogen meter-3'},
            'zooplankton':   {'long_name': 'zooplankton biomass',    'units': 'millimole nitrogen meter-3'},
            'detritus':      {'long_name': 'mole_concentration_of_detritus_expressed_as_nitrogen_in_sea_water',
                              'units': 'millimole nitrogen meter-3'},
        },
        'Fennel': {
            'NO3':      {'long_name': 'nitrate concentration', 'units': 'millimole nitrogen meter-3'},
            'NH4':      {'long_name': 'ammonium concentration', 'units': 'millimole nitrogen meter-3'},
            'PO4':      {'long_name': 'ammonium concentration', 'units': 'millimole po4 meter-3'},
            'chlorophyll': {'long_name': 'chlorophyll concentration', 'units': 'millimole chlorphyll meter-3'},
            'phytoplankton': {'long_name': 'small phytoplankton biomass', 'units': 'millimole nitrogen meter-3'},
            'zooplankton':   {'long_name': 'zooplankton biomass', 'units': 'millimole nitrogen meter-3'},
            'oxygen':    {'long_name': 'oxygen concentration', 'units': 'millimole oxygen meter-3'},
            'TIC':       {'long_name': 'total inorganic carbon', 'units': 'millimole carbon meter-3'},
            'alkalinity':{'long_name': 'total alkalinity', 'units': 'milliequivalent meter-3'},
            'SdetritusC':{'long_name': 'small carbon-detritus concentration', 'units': 'millimole carbon meter-3'},
            'LdetritusC':{'long_name': 'large carbon-detritus concentration', 'units': 'millimole carbon meter-3'},
            'RdetritusC':{'long_name': 'river carbon-detritus concentration', 'units': 'millimole carbon meter-3'},
            'SdetritusN':{'long_name': 'small nitrogen-detritus concentration', 'units': 'millimole nitrogen meter-3'},
            'LdetritusN':{'long_name': 'large nitrogen-detritus concentration', 'units': 'millimole nitrogen meter-3'},
            'RdetritusN':{'long_name': 'river nitrogen-detritus concentration', 'units': 'millimole nitrogen meter-3'}
        }
    }


    if bio_model in bio_tracers:
        for name, attrs in bio_tracers.get(bio_model, {}).items():
            base_variables[name] = ('f4', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), attrs)
    else:
        print("=== Deactivate bio variables ===")


    # NetCDF 생성
    mode = 'w' if cfg.force_write else 'x'
    try:
        ncfile = Dataset(cfg.ininame, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"[✗] {cfg.ininame} already exists and force_write=False")
        return 1
    

    for dim_name, dim_size in dimensions.items():
        ncfile.createDimension(dim_name, dim_size)

    for var_name, (dtype, dims, attrs) in base_variables.items():
        var = ncfile.createVariable(var_name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)

    ncfile.title = cfg.global_attrs.title
    ncfile.clim_file = cfg.ininame
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs.type
    ncfile.history = cfg.global_attrs.type

    sc_r, Cs_r = ru.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ru.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    # 값 할당
    ncfile['spherical'][:] = 'T'
    ncfile['Vtransform'][:] = vtransform
    ncfile['Vstretching'][:] = vstretching
    ncfile['theta_s'][:] = theta_s
    ncfile['theta_b'][:] = theta_b
    ncfile['Tcline'][:] = tcline
    ncfile['hc'][:] = tcline
    ncfile['sc_r'][:] = sc_r
    ncfile['Cs_r'][:] = Cs_r
    ncfile['sc_w'][:] = sc_w
    ncfile['Cs_w'][:] = Cs_w
    ncfile['ocean_time'][0] = initime_num

    for var in ['u', 'v', 'ubar', 'vbar', 'zeta', 'temp', 'salt'] + list(bio_tracers.get(bio_model, {}).keys()):
        ncfile[var][0] = 0

    ncfile.close()

    return 0







