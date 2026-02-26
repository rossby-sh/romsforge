
from __future__ import annotations
from pathlib import Path
import datetime as dt
import numpy as np
import xarray as xr

import yaml


# --------------------------------------------------
# config 자동 로드 (이 스크립트와 같은 디렉토리)
# --------------------------------------------------
HERE = Path(__file__).resolve().parent
CFG_PATH = HERE / "config.yaml"

if not CFG_PATH.exists():
    raise SystemExit(f"config.yaml not found: {CFG_PATH}")

with CFG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)


# --------------------------------------------------
# 유틸
# --------------------------------------------------
def week_start_date_str(t, week_start: str) -> str:
    # weekday(): Monday=0 ... Sunday=6
    wd = t.weekday()
    if week_start == "SUN":
        wd = (wd + 1) % 7
    return (t - dt.timedelta(days=int(wd))).strftime("%Y-%m-%d")


def fmt_monthly_name(template: str, year: int, month: int) -> str:
    return (template
            .replace("{YYYY}", f"{year:04d}")
            .replace("{MM}", f"{month:02d}"))


# --------------------------------------------------
# 경로 / 설정
# --------------------------------------------------
paths = cfg["paths"]
time_cfg = cfg["time"]
opt = cfg.get("options", {})

indir = (HERE / paths.get("indir", ".")).resolve()
outdir = (HERE / paths.get("outdir", paths.get("indir", "."))).resolve()
outdir.mkdir(parents=True, exist_ok=True)

file_list = paths.get("file_list", "avg_use_list.txt")
flist = [
    str(indir / line.strip())
    for line in (indir / file_list).read_text().splitlines()
    if line.strip()
]

year = int(time_cfg["year"])
month = int(time_cfg["month"])
week_start = time_cfg.get("week_start", "MON").upper()
if week_start not in ("MON", "SUN"):
    raise SystemExit("week_start must be MON or SUN")

drop_dup = bool(opt.get("drop_duplicate_time", True))
decode_times = bool(opt.get("decode_times", True))
keep_attrs = bool(opt.get("keep_attrs", True))
write_weekly = bool(opt.get("write_weekly", True))
write_monthly = bool(opt.get("write_monthly", True))
weekly_prefix = opt.get("weekly_prefix", "week")
monthly_template = opt.get("monthly_name", "monthly_{YYYY}{MM}.nc")


# --------------------------------------------------
# 데이터 로드
# --------------------------------------------------
ds = xr.open_mfdataset(
    flist,
    combine="nested",
    concat_dim="ocean_time",
    decode_times=decode_times,
    parallel=False,
).sortby("ocean_time")

if drop_dup:
    tvals = ds["ocean_time"].values
    _, idx = np.unique(tvals, return_index=True)
    ds = ds.isel(ocean_time=np.sort(idx))

time = ds["ocean_time"]

# --------------------------------------------------
# 대상 월 선택 (달력 기준)
# --------------------------------------------------
sel = ds.where(
    (time.dt.year == year) & (time.dt.month == month),
    drop=True,
)

if sel.sizes.get("ocean_time", 0) == 0:
    raise SystemExit(f"No data for {year:04d}-{month:02d}")


# --------------------------------------------------
# week_start 라벨 생성
# --------------------------------------------------
vec = np.vectorize(lambda tt: week_start_date_str(tt, week_start))
wk_label = vec(sel["ocean_time"].values)
sel = sel.assign_coords(week_start=("ocean_time", wk_label))


# --------------------------------------------------
# 주간 평균
# --------------------------------------------------
if write_weekly:
    g = sel.groupby("week_start").mean(dim="ocean_time", keep_attrs=keep_attrs)
    wk_keys = sorted(g["week_start"].values.tolist())

    for i, wk in enumerate(wk_keys, start=1):
        wds = g.sel(week_start=wk)

        src = sel.where(sel["week_start"] == wk, drop=True)
        mid = src["ocean_time"].isel(ocean_time=src.sizes["ocean_time"] // 2)
        wds = wds.expand_dims({"ocean_time": [mid.values]})

        out = outdir / f"{weekly_prefix}{i:02d}.nc"
        wds.to_netcdf(out)
        print("Wrote", out, "(week_start:", wk, ")")


# --------------------------------------------------
# 월평균
# --------------------------------------------------
if write_monthly:
    mmean = sel.mean(dim="ocean_time", keep_attrs=keep_attrs)
    mid = sel["ocean_time"].isel(ocean_time=sel.sizes["ocean_time"] // 2)
    mmean = mmean.expand_dims({"ocean_time": [mid.values]})

    outm = outdir / fmt_monthly_name(monthly_template, year, month)
    mmean.to_netcdf(outm)
    print("Wrote", outm)
