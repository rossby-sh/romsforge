import numpy as np
import datetime as dt
from netCDF4 import Dataset, date2num
import yaml


class ConfigObject:
    def __init__(self, **entries):
        for key, value in entries.items():
            if isinstance(value, dict):
                value = ConfigObject(**value)
            self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self, indent=0):
        lines = []
        pad = '  ' * indent
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigObject):
                lines.append(f"{pad}{key}:")
                lines.append(value.__repr__(indent + 1))
            else:
                lines.append(f"{pad}{key}: {value}")
        return '\n'.join(lines)

    def to_dict(self):
        out = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigObject):
                out[key] = value.to_dict()
            else:
                out[key] = value
        return out

def parse_config(path="config.yaml") -> ConfigObject:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return ConfigObject(**raw)

def compute_relative_time(time, input_unit, ref_unit):
    """
    시간 데이터를 ref_unit 기준 상대 시간으로 변환

    Parameters:
        time        : 원본 시간 배열 (np.ndarray)
        input_unit  : OGCM 시간 단위 (ex: 'hours since ...', 'seconds since ...')
        ref_unit    : ROMS 기준 시간 단위 (ex: 'days since ...')

    Returns:
        np.ndarray : ref_unit 기준 상대 시간
    """
    anchor = dt.datetime(2025, 1, 1)

    # 기준 시간 계산
    t_in  = date2num(anchor, input_unit)
    t_ref = date2num(anchor, ref_unit)

    # 둘 다 days 기준으로 변환
    def to_days(value, unit):
        if "seconds" in unit:
            return value / 86400.0
        elif "hours" in unit:
            return value / 24.0
        elif "days" in unit:
            return value
        else:
            raise ValueError(f"Unknown time unit: {unit}")

    t_in_days = to_days(t_in, input_unit)
    t_ref_days = to_days(t_ref, ref_unit)
    offset = t_in_days - t_ref_days

    # time도 days로 변환
    time_in_days = to_days(time, input_unit)
    total_days = time_in_days + offset

    # 출력 단위 맞추기
    if "seconds" in ref_unit:
        return total_days * 86400.0
    elif "hours" in ref_unit:
        return total_days * 24.0
    elif "days" in ref_unit:
        return total_days
    else:
        raise ValueError(f"Unknown ref_unit: {ref_unit}")



def crop_to_model_domain(lat_src, lon_src, lat_dst, lon_dst):
    """
    lat/lon이 1D (rectilinear) 또는 2D (curvilinear, WRF)인 경우 모두 지원.

    Returns:
        idx, idy: x, y 방향 인덱스 범위
    """

    # 1D version
    if lat_src.ndim == 1 and lon_src.ndim == 1:
        lon_coarse = np.where((lon_src >= np.min(lon_dst)) & (lon_src <= np.max(lon_dst)))[0]
        lat_coarse = np.where((lat_src >= np.min(lat_dst)) & (lat_src <= np.max(lat_dst)))[0]

        lon_ref = lon_src[lon_coarse]
        lat_ref = lat_src[lat_coarse]

        lon_res = max(
            np.max(np.abs(np.diff(lon_ref))),
            np.max(np.abs(np.diff(lon_dst[0, :])))
        )
        lat_res = max(
            np.max(np.abs(np.diff(lat_ref))),
            np.max(np.abs(np.diff(lat_dst[:, 0])))
        )

        lon_min = min(np.min(lon_ref), np.min(lon_dst) - lon_res)
        lon_max = max(np.max(lon_ref), np.max(lon_dst) + lon_res)
        lat_min = min(np.min(lat_ref), np.min(lat_dst) - lat_res)
        lat_max = max(np.max(lat_ref), np.max(lat_dst) + lat_res)

        idx = np.where((lon_src >= lon_min) & (lon_src <= lon_max))[0]
        idy = np.where((lat_src >= lat_min) & (lat_src <= lat_max))[0]

        return idx, idy

    elif lat_src.ndim == 2 and lon_src.ndim == 2:
        lat_min, lat_max = lat_dst.min(), lat_dst.max()
        lon_min, lon_max = lon_dst.min(), lon_dst.max()

        # 해상도 계산 (방향 유의)
        lat_res = np.max(np.abs(np.diff(lat_src[:, 0])))
        lon_res = np.max(np.abs(np.diff(lon_src[0, :])))

        # 확장 경계 계산
        lat_min -= lat_res
        lat_max += lat_res
        lon_min -= lon_res
        lon_max += lon_res

        mask = (
            (lat_src >= lat_min) & (lat_src <= lat_max) &
            (lon_src >= lon_min) & (lon_src <= lon_max)
        )

        j_idx, i_idx = np.where(mask)
        j_min, j_max = j_idx.min(), j_idx.max() + 1
        i_min, i_max = i_idx.min(), i_idx.max() + 1

        idy = np.arange(j_min, j_max)
        idx = np.arange(i_min, i_max)

        return idx, idy


def load_roms_grid(grdname):
    """
    ROMS 그리드 파일에서 필요한 변수들만 추출해 dict 형태로 반환

    Returns:
        dict: {"lon", "lat", "angle", "h", "mask"}
    """
    with Dataset(grdname) as nc:
        lon = nc["lon_rho"][:]
        lat = nc["lat_rho"][:]
        angle = nc["angle"][:]
        h = nc["h"][:]
        mask = nc["mask_rho"][:]
    return ConfigObject(
        lon     = lon,
        lat     = lat,
        angle   = angle,
        topo    = h,
        mask    = mask
    )

def load_ogcm_metadata(ogcm_name,ogcm_var_name):
    """
    OGCM 파일에서 위경도, 깊이, 시간 정보를 읽어와 dict로 반환

    Parameters:
        ogcm_path (str): OGCM NetCDF 파일 경로
        varmap (dict): config.yaml에서 정의한 변수명 매핑

    Returns:
        dict: {"lon", "lat", "depth", "time"}
    """
    with Dataset(ogcm_name) as nc:
        lon = nc[ogcm_var_name["longitude"]][:]
        lat = nc[ogcm_var_name["latitude"]][:]
        depth = nc[ogcm_var_name["depth"]][:]
        time = nc[ogcm_var_name["time"]][:]
        time_unit = nc[ogcm_var_name["time"]].units
    return ConfigObject(
        lon       = lon,
        lat       = lat,
        depth     = depth,
        time      = time,
        time_unit = time_unit
    )

















