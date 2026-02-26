
from pathlib import Path
import yaml
import shutil

HERE = Path(__file__).parent
cfg = yaml.safe_load((HERE / "config_proc.yaml").read_text())

# config 기준으로 매핑
src = Path(cfg["select"]["src_dir"])
dst = Path(cfg["paths"]["outdir"])
dst.mkdir(parents=True, exist_ok=True)

prefix    = cfg["select"]["prefix"]
pattern   = cfg["select"].get("pattern", "{prefix}_avg_{cycle}_{num}.nc")
method    = cfg["select"].get("method", "symlink")
first_use = cfg["select"]["first_use"]
next_use  = cfg["select"]["next_use"]

# cycle 목록 추출: prefix_avg_XXXX_YYYY.nc 에서 XXXX만 뽑는다는 기존 로직 유지
cycles = sorted({p.name.split("_")[2] for p in src.glob(f"{prefix}_avg_*.nc")})

use = []
for i, c in enumerate(cycles):
    nums = first_use if i == 0 else next_use
    for n in nums:
        name = pattern.format(prefix=prefix, cycle=c, num=n)
        f = src / name
        if f.exists():
            use.append(f)

# 기록 + (복사 or 심링크)
list_path = dst / cfg["paths"].get("file_list", "avg_use_list.txt")
with list_path.open("w") as fw:
    for f in use:
        fw.write(f.name + "\n")
        out = dst / f.name

        # 이미 있으면 덮어쓰기 위해 제거
        if out.exists() or out.is_symlink():
            out.unlink()

        if method == "copy":
            shutil.copy2(f, out)
        elif method == "symlink":
            out.symlink_to(f)
        else:
            raise ValueError(f"Unknown method: {method}")

print(f"Wrote list: {list_path}")
print(f"Selected: {len(use)} files")
