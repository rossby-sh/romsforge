# create.py
# - bio 변수 정의를 YAML(bio_vars.yml)에서 읽어서 ini/bry에 각각 다르게 적용
# - base/standard 물리 변수/차원/방향 확장 로직은 기존 그대로 유지
import sys
import os
import numpy as np
import datetime as dt
import yaml
from netCDF4 import Dataset, num2date, date2num
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'libs')))
import datetime as dt
import utils as ut

def load_bio_yaml(path: str) -> dict:
    try:
        with open(path, "r") as f:
            return yaml.safe_load(f) or {}
    except FileNotFoundError:
        # bio yaml이 없으면 bio는 그냥 비활성으로 간주
        return {}

def get_bio_defs(bio_db: dict, bio_model, target: str) -> dict:
    """
    target: 'ini' or 'bry'
    return: {varname: {long_name:..., units:...}, ...}
    """
    if bio_model is None:
        return {}
    key = str(bio_model).lower()
    models = (bio_db or {}).get("bio_models", {})
    if key not in models:
        return {}
    defs = models[key].get(target, {})
    return defs or {}


def create_ini(cfg, grd, initime_num, bio_model=None, ncFormat="NETCDF3_CLASSIC", bio_yaml="bio_vars.yml"):
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask_rho == 1])
    if vtransform == 1 and tcline > hmin_:
        print(f"--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    dimensions = {
        "xi_u": L,
        "xi_v": Lp,
        "xi_rho": Lp,
        "eta_u": Mp,
        "eta_v": M,
        "eta_rho": Mp,
        "s_rho": N,
        "s_w": Np,
        "one": 1,
        "ocean_time": None,
    }

    base_variables = {
        "spherical": ("S1", ("one",), {"long_name": "spherical grid flag"}),
        "Vtransform": ("f4", ("one",), {"long_name": "vertical terrain-following transformation equation"}),
        "Vstretching": ("f4", ("one",), {"long_name": "vertical terrain-following stretching function"}),
        "theta_s": ("f4", ("one",), {"long_name": "S-coordinate surface control parameter", "units": "nondimensional"},),
        "theta_b": ("f4", ("one",), {"long_name": "S-coordinate bottom control parameter", "units": "nondimensional"},),
        "Tcline": ("f4", ("one",), {"long_name": "S-coordinate surface/bottom layer width", "units": "meter"}),
        "hc": ("f4", ("one",), {"long_name": "S-coordinate parameter, critical depth", "units": "meter"}),
        "sc_r": ("f4", ("s_rho",), {"long_name": "S-coordinate at RHO-points", "units": "nondimensional"}),
        "Cs_r": ("f4", ("s_rho",), {"long_name": "S-coordinate stretching curves at RHO-points", "units": "nondimensional"},),
        "sc_w": ("f4", ("s_w",), {"long_name": "S-coordinate at W-points", "units": "nondimensional"}),
        "Cs_w": ("f4",("s_w",),{"long_name": "S-coordinate stretching curves at W-points", "units": "nondimensional"},),
        "ocean_time": ("f4", ("ocean_time",), {"units": cfg.time_ref}),
    }

    init_variables = {
        "u": ("f4", ("ocean_time", "s_rho", "eta_u", "xi_u"), {"long_name": "u-momentum", "units": "meter second-1"}),
        "v": ("f4", ("ocean_time", "s_rho", "eta_v", "xi_v"), {"long_name": "v-momentum", "units": "meter second-1"}),
        "ubar": ("f4", ("ocean_time", "eta_u", "xi_u"), {"long_name": "ubar", "units": "meter second-1"}),
        "vbar": ("f4", ("ocean_time", "eta_v", "xi_v"), {"long_name": "vbar", "units": "meter second-1"}),
        "zeta": ("f4", ("ocean_time", "eta_rho", "xi_rho"), {"long_name": "free surface", "units": "meter"}),
        "temp": ("f4", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), {"long_name": "temperature", "units": "Celsius"},),
        "salt": ("f4", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), {"long_name": "salinity", "units": "PSU"},),
    }

    # --- BIO: YAML에서 ini용 변수 정의 로드 ---
    bio_db = load_bio_yaml(bio_yaml)
    bio_ini = get_bio_defs(bio_db, bio_model, target="ini")

    if bio_ini:
        for name, attrs in bio_ini.items():
            init_variables[name] = ("f4", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), attrs)
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")
    else:
        print("--- [NOTE] Deactivate initiating biological variables ---")

    mode = "w" if cfg.force_write else "x"
    try:
        ncfile = Dataset(cfg.ininame, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {cfg.ininame} already exists and force_write=False ---")
        return 1

    for dim_name, dim_size in dimensions.items():
        ncfile.createDimension(dim_name, dim_size)

    for name, (dtype, dims, attrs) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)

    # init_variables는 껍데기 만들 때 0으로 초기화(너 정책 유지)
    for name, (dtype, dims, attrs) in init_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)
        var[0] = 0.0

    sc_r, Cs_r = ut.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ut.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    ncfile["sc_r"][:] = sc_r
    ncfile["Cs_r"][:] = Cs_r
    ncfile["sc_w"][:] = sc_w
    ncfile["Cs_w"][:] = Cs_w
    ncfile["theta_s"][:] = theta_s
    ncfile["theta_b"][:] = theta_b
    ncfile["Tcline"][:] = tcline
    ncfile["hc"][:] = tcline
    ncfile["Vtransform"][:] = vtransform
    ncfile["Vstretching"][:] = vstretching
    ncfile["spherical"][:] = "T"
    ncfile["ocean_time"][0] = initime_num

    ncfile.title = cfg.global_attrs.title
    ncfile.clim_file = cfg.ininame
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs.type
    ncfile.history = f"Created on {dt.datetime.now().isoformat()}"

    ncfile.close()
    print(f"--- [+] Initial file created: {cfg.ininame} ---")
    return 0


def create_bry(cfg, grd, bry_time, bio_model=None, ncFormat="NETCDF3_CLASSIC", bio_yaml="bio_vars.yml"):
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask_rho == 1])
    if vtransform == 1 and tcline > hmin_:
        print("--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    directions = ["north", "south", "east", "west"]
    grids = {"north": "xi", "south": "xi", "east": "eta", "west": "eta"}  # (현재 코드 유지용, 안 쓰면 지워도 됨)

    dimensions = {
        "xi_u": L,
        "xi_v": Lp,
        "xi_rho": Lp,
        "eta_u": Mp,
        "eta_v": M,
        "eta_rho": Mp,
        "s_rho": N,
        "s_w": Np,
        "one": 1,
        "bry_time": len(bry_time),
        "zeta_time": len(bry_time),
        "temp_time": len(bry_time),
        "salt_time": len(bry_time),
        "v2d_time": len(bry_time),
        "v3d_time": len(bry_time),
    }

    # --- base_variables: stretching 포함 (기존 그대로) ---
    sc_r, Cs_r = ut.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ut.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    base_variables = {
        "spherical": ("S1", ("one",), {}, "T"),
        "Vtransform": (
            "f4",
            ("one",),
            {"long_name": "vertical terrain-following transformation equation"},
            vtransform,
        ),
        "Vstretching": (
            "f4",
            ("one",),
            {"long_name": "vertical terrain-following stretching function"},
            vstretching,
        ),
        "theta_s": (
            "f4",
            ("one",),
            {"long_name": "S-coordinate surface control parameter", "units": "nondimensional"},
            theta_s,
        ),
        "theta_b": (
            "f4",
            ("one",),
            {"long_name": "S-coordinate bottom control parameter", "units": "nondimensional"},
            theta_b,
        ),
        "Tcline": (
            "f4",
            ("one",),
            {"long_name": "S-coordinate surface/bottom layer width", "units": "meter"},
            tcline,
        ),
        "hc": (
            "f4",
            ("one",),
            {"long_name": "S-coordinate parameter, critical depth", "units": "meter"},
            tcline,
        ),
        "sc_r": ("f4", ("s_rho",), {"long_name": "S-coordinate at RHO-points", "units": "nondimensional"}, sc_r),
        "Cs_r": (
            "f4",
            ("s_rho",),
            {"long_name": "S-coordinate stretching curves at RHO-points", "units": "nondimensional"},
            Cs_r,
        ),
        "sc_w": ("f4", ("s_w",), {"long_name": "S-coordinate at W-points", "units": "nondimensional"}, sc_w),
        "Cs_w": (
            "f4",
            ("s_w",),
            {"long_name": "S-coordinate stretching curves at W-points", "units": "nondimensional"},
            Cs_w,
        ),
    }

    time_dims = ["bry_time", "zeta_time", "temp_time", "salt_time", "v2d_time", "v3d_time"]

    standard_vars = {
        "zeta": ("f4", ("zeta_time", "eta_rho", "xi_rho"), {"long_name": "free-surface", "units": "meter"}),
        "ubar": (
            "f4",
            ("v2d_time", "eta_u", "xi_u"),
            {"long_name": "vertically integrated u-momentum", "units": "meter second-1"},
        ),
        "vbar": (
            "f4",
            ("v2d_time", "eta_v", "xi_v"),
            {"long_name": "vertically integrated v-momentum", "units": "meter second-1"},
        ),
        "u": ("f4", ("v3d_time", "s_rho", "eta_u", "xi_u"), {"long_name": "u-momentum", "units": "meter second-1"}),
        "v": ("f4", ("v3d_time", "s_rho", "eta_v", "xi_v"), {"long_name": "v-momentum", "units": "meter second-1"}),
        "temp": (
            "f4",
            ("temp_time", "s_rho", "eta_rho", "xi_rho"),
            {"long_name": "potential temperature", "units": "Celsius"},
        ),
        "salt": ("f4", ("salt_time", "s_rho", "eta_rho", "xi_rho"), {"long_name": "salinity", "units": "PSU"}),
    }

    # --- BIO: YAML에서 bry용 변수 정의 로드 ---
    bio_db = load_bio_yaml(bio_yaml)
    bio_bry = get_bio_defs(bio_db, bio_model, target="bry")

    all_vars = dict(standard_vars)
    if bio_bry:
        # bry는 아래에서 방향별 확장을 하기 때문에, 일단 full dims로 넣고 get_bry_dims가 잘라먹게 둔다.
        for name, attrs in bio_bry.items():
            all_vars[name] = ("f4", ("bry_time", "s_rho", "eta_rho", "xi_rho"), attrs)
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")
    else:
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type ---")

    mode = "w" if cfg.force_write else "x"
    try:
        ncfile = Dataset(cfg.bryname, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {cfg.bryname} already exists and force_write=False ---")
        return 1

    for name, size in dimensions.items():
        ncfile.createDimension(name, size)

    for name, (dtype, dims, attrs, value) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = value

    def get_bry_dims(dims, direction):
        keep_dim = "xi_" if direction in ["north", "south"] else "eta_"
        return tuple(d for d in dims if not d.startswith(("xi_", "eta_")) or d.startswith(keep_dim))

    bry_variables = {}
    for varname, (dtype, dims, attrs) in all_vars.items():
        for d in directions:
            dims_d = get_bry_dims(dims, d)
            bry_variables[f"{varname}_{d}"] = (
                dtype,
                dims_d,
                {**attrs, "long_name": f"{attrs['long_name']} at {d} boundary"},
            )

    # Time variables
    for name in time_dims:
        v = ncfile.createVariable(name, "f4", (name,))
        v.units = cfg.time_ref
        v.long_name = f"time for {name.replace('_time','')} condition"
        v[:] = bry_time

    # Boundary variables (껍데기니까 0으로 채움)
    for name, (dtype, dims, attrs) in bry_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        var.setncatts(attrs)
        var[:] = 0.0

    ncfile.title = cfg.global_attrs_bry.title
    ncfile.clim_file = cfg.bryname
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs_bry.type
    ncfile.history = f"Created on {dt.datetime.now().isoformat()}"

    ncfile.close()
    print(f"--- [+] Boundary file created: {cfg.bryname} ---")
    return 0


def create_std(cfg, grd, stdtime_num, outname, bio_model=None, ncFormat="NETCDF3_CLASSIC", bio_yaml="bio_vars.yml"):
    # NOTE: 구조/스타일은 create_ini 그대로 유지하고, 출력 파일만 outname으로 분리
    vstretching, vtransform = cfg.vertical.vstretching, cfg.vertical.vtransform
    theta_s, theta_b = cfg.vertical.theta_s, cfg.vertical.theta_b
    tcline, layer_n = cfg.vertical.tcline, cfg.vertical.layer_n

    hmin_ = np.min(grd.topo[grd.mask_rho == 1])
    if vtransform == 1 and tcline > hmin_:
        print(f"--- [!ERROR] Tcline must be <= hmin when Vtransform == 1 ---")
        return 1

    Mp, Lp = grd.topo.shape
    L, M, N, Np = Lp - 1, Mp - 1, layer_n, layer_n + 1

    dimensions = {
        "xi_u": L,
        "xi_v": Lp,
        "xi_rho": Lp,
        "eta_u": Mp,
        "eta_v": M,
        "eta_rho": Mp,
        "s_rho": N,
        "s_w": Np,
        "one": 1,
        "ocean_time": None,
    }

    base_variables = {
        "spherical": ("S1", ("one",), {"long_name": "spherical grid flag"}),
        "Vtransform": ("f4", ("one",), {"long_name": "vertical terrain-following transformation equation"}),
        "Vstretching": ("f4", ("one",), {"long_name": "vertical terrain-following stretching function"}),
        "theta_s": ("f4", ("one",), {"long_name": "S-coordinate surface control parameter", "units": "nondimensional"},),
        "theta_b": ("f4", ("one",), {"long_name": "S-coordinate bottom control parameter", "units": "nondimensional"},),
        "Tcline": ("f4", ("one",), {"long_name": "S-coordinate surface/bottom layer width", "units": "meter"}),
        "hc": ("f4", ("one",), {"long_name": "S-coordinate parameter, critical depth", "units": "meter"}),
        "sc_r": ("f4", ("s_rho",), {"long_name": "S-coordinate at RHO-points", "units": "nondimensional"}),
        "Cs_r": ("f4", ("s_rho",), {"long_name": "S-coordinate stretching curves at RHO-points", "units": "nondimensional"},),
        "sc_w": ("f4", ("s_w",), {"long_name": "S-coordinate at W-points", "units": "nondimensional"}),
        "Cs_w": ("f4",("s_w",),{"long_name": "S-coordinate stretching curves at W-points", "units": "nondimensional"},),
        "ocean_time": ("f4", ("ocean_time",), {"units": cfg.time_ref}),
    }

    std_variables = {
        "u": ("f4", ("ocean_time", "s_rho", "eta_u", "xi_u"), {"long_name": "u-momentum", "units": "meter second-1"}),
        "v": ("f4", ("ocean_time", "s_rho", "eta_v", "xi_v"), {"long_name": "v-momentum", "units": "meter second-1"}),
        "ubar": ("f4", ("ocean_time", "eta_u", "xi_u"), {"long_name": "ubar", "units": "meter second-1"}),
        "vbar": ("f4", ("ocean_time", "eta_v", "xi_v"), {"long_name": "vbar", "units": "meter second-1"}),
        "zeta": ("f4", ("ocean_time", "eta_rho", "xi_rho"), {"long_name": "free surface", "units": "meter"}),
        "temp": ("f4", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), {"long_name": "temperature", "units": "Celsius"},),
        "salt": ("f4", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), {"long_name": "salinity", "units": "PSU"},),
    }

    # --- [bio variables] (create_ini와 동일: YAML에서 bio 변수 정의 로드) ---
    # 1) bio_model 우선순위: 인자 > cfg.bio_type > cfg.bio_model_type

    if bio_model is None:
        bio_model = getattr(cfg, "bio_model_type", None)

    bio_db = load_bio_yaml(bio_yaml)

    # 2) std 타겟 정의가 있으면 그걸 쓰고, 없으면 ini 정의를 fallback으로 사용
    bio_std = get_bio_defs(bio_db, bio_model, target="std")
    if not bio_std:
        bio_std = get_bio_defs(bio_db, bio_model, target="ini")
        if bio_std:
            print(f"--- [NOTE] Use 'ini' biological variable definitions for std target: {bio_model} ---")

    if bio_std:
        for name, attrs in bio_std.items():
            std_variables[name] = ("f4", ("ocean_time", "s_rho", "eta_rho", "xi_rho"), attrs)
        print(f"--- [NOTE] Initiating biological variables: {bio_model} type (std) ---")
    else:
        if bio_model is not None:
            print(f"--- [NOTE] No biological variables found for std target in YAML: {bio_model} ---")
        else:
            print("--- [NOTE] Deactivate initiating biological variables ---")

    outdir = os.path.dirname(outname)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    mode = "w" if cfg.force_write else "x"
    try:
        ncfile = Dataset(outname, mode=mode, format=ncFormat)
    except FileExistsError:
        print(f"--- [!ERROR] {outname} already exists and force_write=False ---")
        return 1

    for dim_name, dim_size in dimensions.items():
        ncfile.createDimension(dim_name, dim_size)

    for name, (dtype, dims, attrs) in base_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)

    # std_variables는 껍데기 만들 때 0으로 초기화(ini 정책 유지)
    for name, (dtype, dims, attrs) in std_variables.items():
        var = ncfile.createVariable(name, dtype, dims)
        for attr_name, attr_val in attrs.items():
            setattr(var, attr_name, attr_val)
        var[0] = 0.0

    sc_r, Cs_r = ut.stretching(vstretching, theta_s, theta_b, layer_n, 0)
    sc_w, Cs_w = ut.stretching(vstretching, theta_s, theta_b, layer_n, 1)

    ncfile["sc_r"][:] = sc_r
    ncfile["Cs_r"][:] = Cs_r
    ncfile["sc_w"][:] = sc_w
    ncfile["Cs_w"][:] = Cs_w

    ncfile["hc"][:] = tcline
    ncfile["Tcline"][:] = tcline
    ncfile["theta_s"][:] = theta_s
    ncfile["theta_b"][:] = theta_b
    ncfile["Vtransform"][:] = vtransform
    ncfile["Vstretching"][:] = vstretching
    ncfile["spherical"][:] = "T"
    ncfile["ocean_time"][0] = stdtime_num

    ncfile.title = cfg.global_attrs.title
    ncfile.clim_file = outname
    ncfile.grd_file = cfg.grdname
    ncfile.type = cfg.global_attrs.type
    ncfile.history = f"Created on {dt.datetime.now().isoformat()}"
    ncfile.applied_scale_factor = str(cfg.scale_factor)
    if bio_model is not None:
        ncfile.applied_bio_scale_factor = str(cfg.bio_scale_factor)
    ncfile.std_ogcm_name = cfg.ogcm_name
    ncfile.std_bio_ogcm_name = cfg.bio_ogcm_name
    ncfile.close()
    print(f"--- [+] STD file created: {outname} ---")
    return 0


def create_obs(cfg, grd, outname, survey_time, Nsurvey, ncFormat="NETCDF3_CLASSIC", unlimited_datum=True):
    """
    - survey dimension: fixed (Nsurvey)
    - datum dimension: unlimited if unlimited_datum else fixed (not used now)
    """
    outdir = os.path.dirname(outname)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    mode = "w" if getattr(cfg, "force_write", True) else "x"
    try:
        nc = Dataset(outname, mode=mode, format=("NETCDF4" if "NETCDF4" in ncFormat else ncFormat))
    except FileExistsError:
        return 1

    # --- dims ---
    nc.createDimension("one", 1)
    nc.createDimension("state_variable", int(cfg.Nstate))
    nc.createDimension("survey", int(Nsurvey))
    nc.createDimension("datum", None if unlimited_datum else 0)

    # --- vars ---
    zlib = ("NETCDF4" in nc.file_format)

    v = nc.createVariable("spherical", "S1", ("one",), zlib=zlib)
    v.long_name = "grid type logical switch"
    v.option_T = "spherical"
    v.option_F = "Cartesian"
    v[:] = np.array(["T"], dtype="S1")

    def _v(name, dtype, dims, attrs):
        vv = nc.createVariable(name, dtype, dims, zlib=zlib)
        for k, val in attrs.items():
            setattr(vv, k, val)
        return vv

    _v("obs_type", "i4", ("datum",), {"long_name": "model state variable associated with observation"})
    _v("obs_time", "f8", ("datum",), {"long_name": "Time of observation", "units": "days", "calendar": "standard"})
    _v("obs_lon",  "f8", ("datum",), {"long_name": "Longitude of observation", "units": "degrees_east"})
    _v("obs_lat",  "f8", ("datum",), {"long_name": "Latitude of observation",  "units": "degrees_north"})
    _v("obs_depth","f8", ("datum",), {"long_name": "Depth of observation", "units": "meter", "minus": "downwards"})
    _v("obs_error","f8", ("datum",), {"long_name": "Observation error covariance"})
    _v("obs_value","f8", ("datum",), {"long_name": "Observation value"})
    _v("obs_Xgrid","f8", ("datum",), {"long_name": "x-grid observation location", "units": "nondimensional"})
    _v("obs_Ygrid","f8", ("datum",), {"long_name": "y-grid observation location", "units": "nondimensional"})
    _v("obs_Zgrid","f8", ("datum",), {"long_name": "z-grid observation location", "units": "nondimensional"})
    _v("obs_provenance","f8", ("datum",), {"long_name": "observation origin"})

    _v("survey_time","f8", ("survey",), {"long_name": "Survey time", "units": "day", "calendar": "standard"})
    _v("Nobs","i4", ("survey",), {"long_name": "number of observations with the same survey time"})
    _v("obs_variance","f8", ("state_variable",), {"long_name": "global time and space observation variance"})

    # fill survey_time now (A)
    nc["survey_time"][:] = np.asarray(survey_time, dtype=np.float64)
    nc["Nobs"][:] = 0
    nc["obs_variance"][:] = 0.0

    # attrs
    nc.title = cfg.global_attrs.title
    nc.grd_file = cfg.grdname
    nc.type = cfg.global_attrs.type
    nc.history = f"Created on {dt.datetime.now().isoformat()}"

    nc.close()
    return 0


def append_obs(outname, chunk, ncFormat="NETCDF3_CLASSIC"):
    """
    Append chunk arrays to datum dimension.
    chunk: dict returned by build_obs_chunk_from_source
    """
    with Dataset(outname, mode="a") as nc:
        n0 = nc.dimensions["datum"].size
        n1 = n0 + int(chunk["Nobs"])

        for key in ["obs_type","obs_time","obs_lon","obs_lat","obs_depth","obs_error","obs_value",
                    "obs_Xgrid","obs_Ygrid","obs_Zgrid","obs_provenance"]:
            nc[key][n0:n1] = chunk[key]


def finalize_obs(outname, survey_time, Nobs, obs_variance, ncFormat="NETCDF3_CLASSIC"):
    with Dataset(outname, mode="a") as nc:
        nc["survey_time"][:] = np.asarray(survey_time, dtype=np.float64)
        nc["Nobs"][:] = np.asarray(Nobs, dtype=np.int32)
        nc["obs_variance"][:] = np.asarray(obs_variance, dtype=np.float64)







