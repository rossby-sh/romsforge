
from pathlib import Path
import yaml

HERE = Path(__file__).parent
cfg = yaml.safe_load((HERE / "config_proc.yaml").read_text())

src = Path(cfg["paths"]["src_dir"])
dst = Path(cfg["paths"]["post_proc"])
dst.mkdir(exist_ok=True)

prefix = cfg["avg"]["prefix"]
first_use = cfg["avg"]["first_use"]
next_use  = cfg["avg"]["next_use"]

# cycle 목록
cycles = sorted({p.name.split("_")[2] for p in src.glob(f"{prefix}_avg_*.nc")})

use = []
for i, c in enumerate(cycles):
    nums = first_use if i == 0 else next_use
    for n in nums:
        f = src / f"{prefix}_avg_{c}_{n}.nc"
        if f.exists():
            use.append(f)

# 기록 + 복사
with (dst / "avg_use_list.txt").open("w") as fw:
    for f in use:
        fw.write(f.name + "\n")
        (dst / f.name).symlink_to(f)
