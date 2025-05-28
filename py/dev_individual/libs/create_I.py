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
import roms_utils as ru  # stretching 함수가 포함된 모듈

def create_ini(ininame, mask, topo, MyVar, initime_num, time_ref, Title,
                    ncFormat='NETCDF3_CLASSIC', bio_model=None):

    hmin_ = np.min(topo[mask == 1])
    if MyVar['Vtransform'] == 1 and MyVar['Tcline'] > hmin_:
        raise ValueError("Tcline must be <= hmin when Vtransform == 1")

    Mp, Lp = topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, MyVar['Layer_N'], MyVar['Layer_N'] + 1

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

        'ocean_time':    ('f4',   ('ocean_time',), {'units': time_ref}),

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
            'phytoplankton': {'long_name': 'small phytoplankton biomass', 'units': 'millimole nitrogen meter-3'},
            'zooplankton':   {'long_name': 'zooplankton biomass', 'units': 'millimole nitrogen meter-3'},
            'Ldetritus': {'long_name': 'large detritus nitrogen concentration', 'units': 'millimole nitrogen meter-3'},
            'Sdetritus': {'long_name': 'small detritus nitrogen concentration', 'units': 'millimole nitrogen meter-3'},
            'oxygen':    {'long_name': 'oxygen concentration', 'units': 'millimole oxygen meter-3'},
            'TIC':       {'long_name': 'total inorganic carbon', 'units': 'millimole carbon meter-3'},
            'alkalinity':{'long_name': 'total alkalinity', 'units': 'milliequivalent meter-3'},
        }
    }


    if bio_model in bio_tracers:
        for name, attrs in bio_tracers.get(bio_model, {}).items():
            base_variables[name] = ('f4', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), attrs)
    else:
        print("=== Deactivate bio variables ===")


    # NetCDF 생성
    ncfile = Dataset(ininame, mode='w', format=ncFormat)
    for dim_name, dim_size in dimensions.items():
        ncfile.createDimension(dim_name, dim_size)

    for var_name, (dtype, dims, attrs) in base_variables.items():
        var = ncfile.createVariable(var_name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)

    ncfile.title = Title
    ncfile.clim_file = ininame
    ncfile.grd_file = ''
    ncfile.type = 'INITIAL file'
    ncfile.history = 'ROMS'

    sc_r, Cs_r = ru.stretching(MyVar['Vstretching'], MyVar['Theta_s'], MyVar['Theta_b'], MyVar['Layer_N'], 0)
    sc_w, Cs_w = ru.stretching(MyVar['Vstretching'], MyVar['Theta_s'], MyVar['Theta_b'], MyVar['Layer_N'], 1)

    # 값 할당
    ncfile['spherical'][:] = 'T'
    ncfile['Vtransform'][:] = MyVar['Vtransform']
    ncfile['Vstretching'][:] = MyVar['Vstretching']
    ncfile['theta_s'][:] = MyVar['Theta_s']
    ncfile['theta_b'][:] = MyVar['Theta_b']
    ncfile['Tcline'][:] = MyVar['Tcline']
    ncfile['hc'][:] = MyVar['Tcline']
    ncfile['sc_r'][:] = sc_r
    ncfile['Cs_r'][:] = Cs_r
    ncfile['sc_w'][:] = sc_w
    ncfile['Cs_w'][:] = Cs_w
    ncfile['ocean_time'][0] = initime_num

    for var in ['u', 'v', 'ubar', 'vbar', 'zeta', 'temp', 'salt'] + list(bio_tracers.get(bio_model, {}).keys()):
        ncfile[var][0] = 0

    ncfile.close()









































