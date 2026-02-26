
#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import datetime as pydt

import numpy as np
import xarray as xr
import yaml
from netCDF4 import num2date, Dataset


# ============================================================
# config auto-load
# ============================================================
HERE = Path(__file__).resolve().parent
CFG_PATH = HERE / "config_proc.yaml"

if not CFG_PATH.exists():
    raise SystemExit(f"config_proc.yaml not found: {CFG_PATH}")

with CFG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)


def fatal(msg: str) -> None:
    raise SystemExit(f"[FATAL] {msg}")


# ============================================================
# selection utils
# ============================================================
def get_cycles(src_dir: Path, prefix: str, kind: str) -> list[str]:
    files = sorted(src_dir.glob(f"{prefix}_{kind}_*.nc"))
    cycles = set()
    for p in files:
        parts = p.name.split("_")
        if len(parts) < 4:
            continue
        cycle = parts[2]
        if cycle.isdigit() and len(cycle) == 4:
            cycles.add(cycle)
    return sorted(cycles, key=lambda x: int(x))


def select_files(cfg: dict) -> tuple[list[Path], str]:
    """
    Returns:
      chosen_files, kind
    """
    sel = cfg.get("select")
    if not sel:
        fatal("Missing 'select' block in config_proc.yaml")

    required = ["src_dir", "prefix", "first_use", "next_use"]
    for k in required:
        if k not in sel:
            fatal(f"Missing select.{k} in config_proc.yaml")

    src_dir = Path(sel["src_dir"]).expanduser().resolve()
    prefix = str(sel["prefix"])

    kind = str(sel.get("kind", "avg")).lower()  # "avg" or "his"
    if kind not in ("avg", "his"):
        fatal("select.kind must be 'avg' or 'his'")

    first_use = sel["first_use"]
    next_use = sel["next_use"]

    # default: {prefix}_{kind}_{cycle}_{num}.nc
    pattern = sel.get("pattern", "{prefix}_{kind}_{cycle}_{num}.nc")

    if not isinstance(first_use, list) or not first_use:
        fatal("select.first_use must be a non-empty list")
    if not isinstance(next_use, list) or not next_use:
        fatal("select.next_use must be a non-empty list")

    if not src_dir.exists():
        fatal(f"select.src_dir does not exist: {src_dir}")

    cycles = get_cycles(src_dir, prefix, kind)
    if not cycles:
        fatal(f"No {kind} files found in {src_dir} with prefix={prefix}")

    chosen: list[Path] = []
    for i, cycle in enumerate(cycles):
        nums = first_use if i == 0 else next_use
        for num in nums:
            name = pattern.format(prefix=prefix, kind=kind, cycle=cycle, num=num)
            p = src_dir / name
            if p.exists():
                chosen.append(p)

    if not chosen:
        fatal("Selection resulted in 0 files. Check select settings.")
    return chosen, kind


def write_file_list(cfg: dict, chosen: list[Path]) -> Path:
    """
    기록용: 선택된 원본 파일 절대경로를 outdir에 저장 (staging 없음).
    """
    paths = cfg.get("paths", {}) or {}
    if "outdir" not in paths:
        fatal("Missing paths.outdir in config_proc.yaml")

    outdir = Path(paths["outdir"]).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    file_list_name = paths.get("file_list", "use_list.txt")
    file_list_path = outdir / file_list_name

    lines = [str(p.resolve()) for p in chosen]
    file_list_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return file_list_path


# ============================================================
# time grouping (calendar week, cftime-safe)
# ============================================================
def to_py_datetime(d) -> pydt.datetime:
    sec = int(getattr(d, "second", 0))
    return pydt.datetime(
        int(d.year), int(d.month), int(d.day),
        int(getattr(d, "hour", 0)),
        int(getattr(d, "minute", 0)),
        sec
    )


def build_week_start_keys(dt_all: np.ndarray, week_start: str) -> np.ndarray:
    py_dt = np.array([to_py_datetime(d) for d in dt_all], dtype=object)

    # python weekday: Monday=0 ... Sunday=6
    wd = np.array([d.weekday() for d in py_dt], dtype=int)

    # SUN start이면 Sunday=0 ... Saturday=6 로 변환
    if week_start == "SUN":
        wd = (wd + 1) % 7

    wk_start_dt = np.array(
        [py_dt[i] - pydt.timedelta(days=int(wd[i])) for i in range(len(py_dt))],
        dtype=object
    )
    return np.array([d.strftime("%Y-%m-%d") for d in wk_start_dt], dtype=object)


def build_target_mask(dt_all: np.ndarray, year: int, month: int) -> np.ndarray:
    py_dt = np.array([to_py_datetime(d) for d in dt_all], dtype=object)
    return np.array([(d.year == year and d.month == month) for d in py_dt], dtype=bool)


# ============================================================
# postprocess - add chlorophyll (netCDF4 append)
# ============================================================
def ensure_var_like(nc: Dataset, src_name: str, dst_name: str, dtype="f4"):
    if src_name not in nc.variables:
        fatal(f"variable '{src_name}' not found in {nc.filepath()}")

    src = nc.variables[src_name]

    if dst_name in nc.variables:
        dst = nc.variables[dst_name]
        return src, dst

    dst = nc.createVariable(dst_name, dtype, src.dimensions, zlib=False)
    dst.long_name = "Chlorophyll-a"
    dst.units = "mg m-3"
    dst.source = "derived from phytoplankton"
    return src, dst


def add_chl_inplace(nc_path: Path, factor: float, overwrite: bool = True):
    with Dataset(nc_path, "a") as nc:
        src_name = "phytoplankton"
        dst_name = "chlorophyll"

        if (dst_name in nc.variables) and (not overwrite):
            print(f"Skip (exists): {nc_path}")
            return

        src, dst = ensure_var_like(nc, src_name, dst_name, dtype="f4")
        dst.factor = float(factor)

        shape = src.shape
        ndim = len(shape)

        if ndim <= 1:
            dst[:] = (np.asarray(src[:], dtype=np.float32) * np.float32(factor))
            return

        # 메모리 보호: axis=1 기준 슬라이스
        chunk_axis = 1 if ndim >= 2 else 0
        n = shape[chunk_axis]
        sl = [slice(None)] * ndim

        for i in range(n):
            sl[chunk_axis] = slice(i, i + 1)
            block = np.asarray(src[tuple(sl)], dtype=np.float32) * np.float32(factor)
            dst[tuple(sl)] = block


# ============================================================
# main
# ============================================================
def main():
    # naming
    naming = cfg.get("naming", {}) or {}
    if "prefix" not in naming:
        fatal("Missing naming.prefix in config_proc.yaml")
    name_prefix = str(naming["prefix"])

    # time
    time_cfg = cfg.get("time", {}) or {}
    for k in ["year", "month", "week_start"]:
        if k not in time_cfg:
            fatal(f"Missing time.{k} in config_proc.yaml")
    year = int(time_cfg["year"])
    month = int(time_cfg["month"])
    week_start = str(time_cfg["week_start"]).upper()
    if week_start not in ("MON", "SUN"):
        fatal("time.week_start must be 'MON' or 'SUN'")

    # options
    opt = cfg.get("options", {}) or {}
    drop_dup = bool(opt.get("drop_duplicate_time", True))
    keep_attrs = bool(opt.get("keep_attrs", True))
    write_weekly = bool(opt.get("write_weekly", True))
    write_monthly = bool(opt.get("write_monthly", True))
    add_chl = bool(opt.get("add_chl", False))

    # paths
    paths = cfg.get("paths", {}) or {}
    if "outdir" not in paths:
        fatal("Missing paths.outdir in config_proc.yaml")
    outdir = Path(paths["outdir"]).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)


    use_existing = bool(opt.get("use_existing_file_list", False))

    paths = cfg.get("paths", {}) or {}
    outdir = Path(paths["outdir"]).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    file_list_name = paths.get("file_list", "use_list.txt")
    file_list_path = outdir / file_list_name

    if use_existing and file_list_path.exists():
        print(f"Using existing file list: {file_list_path}")
        flist = file_list_path.read_text().splitlines()

        if len(flist) == 0:
            fatal(f"File list is empty: {file_list_path}")

        chosen = [Path(p) for p in flist]
        kind = "avg"   # fallback (optional)

    else:
        print("Generating new file list from config selection...")
        chosen, kind = select_files(cfg)
        file_list_path = write_file_list(cfg, chosen)
        flist = [str(p) for p in chosen]

    print(f"Total input files: {len(flist)}")
    print(f"File list path: {file_list_path}")


    ds = xr.open_mfdataset(
        flist,
        combine="nested",
        concat_dim="ocean_time",
        decode_times=False,   # ocean_time 숫자 보존
        parallel=False,
    ).sortby("ocean_time")

    # 3) drop duplicate ocean_time (numeric)
    if drop_dup:
        tvals = ds["ocean_time"].values
        _, idx = np.unique(tvals, return_index=True)
        ds = ds.isel(ocean_time=np.sort(idx))

    # 4) ocean_time meta
    ot = ds["ocean_time"].values
    units = ds["ocean_time"].attrs.get("units", None)
    calendar = ds["ocean_time"].attrs.get("calendar", "standard")
    if units is None:
        fatal("ocean_time:units missing. (check original files)")

    # 5) grouping datetime
    dt_all = np.array(num2date(ot, units=units, calendar=calendar), dtype=object)

    # 6) week key + month mask
    week_start_key = build_week_start_keys(dt_all, week_start)
    target_mask = build_target_mask(dt_all, year, month)
    if not target_mask.any():
        fatal(f"No data in {year:04d}-{month:02d} within selected files ({kind}).")

    # 출력할 주 리스트는 "그 달에 걸치는 주"만
    wk_keys = np.unique(week_start_key[target_mask]).tolist()
    wk_keys = sorted([str(x) for x in wk_keys])
    if len(wk_keys) == 0:
        fatal(f"No weeks intersecting {year:04d}-{month:02d}")

    # groupby coord
    ds = ds.assign_coords(week_start=("ocean_time", week_start_key))

    # encoding
    enc_time = {"ocean_time": {"dtype": "float64", "_FillValue": None}}

    written_files: list[Path] = []

    # ============================================================
    # MONTHLY FIRST (요청: 월평균 먼저 파일 생성)
    # ============================================================
    if write_monthly:
        month_idx = np.where(target_mask)[0]
        if month_idx.size == 0:
            fatal(f"No data inside {year:04d}-{month:02d} for monthly mean")

        sel_month = ds.isel(ocean_time=month_idx)
        mmean = sel_month.mean(dim="ocean_time", keep_attrs=keep_attrs)

        ntime = sel_month.sizes["ocean_time"]
        mid_num = float(sel_month["ocean_time"].values[ntime // 2])
        mmean = mmean.expand_dims({"ocean_time": [mid_num]})

        mmean = mmean.drop_vars("week_start", errors="ignore")
        mmean["ocean_time"].attrs["units"] = units
        mmean["ocean_time"].attrs["calendar"] = calendar

        outm = outdir / f"{name_prefix}_{year:04d}{month:02d}_monthly.nc"
        mmean.to_netcdf(outm, encoding=enc_time)
        written_files.append(outm)
        print(f"Wrote {outm}")

    # ============================================================
    # WEEKLY NEXT (풀 주 로직: 평균은 월 경계 넘어가도 주 전체 사용)
    # ============================================================
    if write_weekly:
        wk_arr = week_start_key
        for i, wk in enumerate(wk_keys, start=1):
            # ★ 풀 주: target_mask로 자르지 않음
            idx = np.where(wk_arr == wk)[0]
            if idx.size == 0:
                fatal(f"Internal error: week_start={wk} has 0 samples")

            src = ds.isel(ocean_time=idx)
            wds = src.mean(dim="ocean_time", keep_attrs=keep_attrs)

            mid_num = float(src["ocean_time"].values[idx.size // 2])
            wds = wds.expand_dims({"ocean_time": [mid_num]})

            wds = wds.drop_vars("week_start", errors="ignore")
            wds["ocean_time"].attrs["units"] = units
            wds["ocean_time"].attrs["calendar"] = calendar

            out = outdir / f"{name_prefix}_{year:04d}{month:02d}_week{i:02d}.nc"
            wds.to_netcdf(out, encoding=enc_time)
            written_files.append(out)
            print(f"Wrote {out}  (week_start={wk}, ntime={idx.size})")

    print(f"Done. kind={kind}")
    print(f"Selected files: {len(chosen)}")
    print(f"File list written: {file_list_path}")
    print(f"Output directory: {outdir}")
    print(f"ocean_time preserved: units='{units}', calendar='{calendar}'")

    # post: add chlorophyll
    if add_chl:
        factor = 0.02 * 6.625 * 12
        print(f"Post: add chlorophyll (factor={factor}) to {len(written_files)} files...")
        for f in written_files:
            add_chl_inplace(f, factor=factor, overwrite=True)
        print("Post: add chlorophyll done.")
    else:
        print("Post: add_chl=false -> skip.")


if __name__ == "__main__":
    main()
