
#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./select_avg_files.sh <SRC_DIR> <DST_DIR> [PREFIX]
#
# Example:
#   ./select_avg_files.sh /path/to/roms/run /path/to/post_proc nifs5km
#
# Output:
#   <DST_DIR>/avg_use_list.txt  (선택된 파일 목록)
#   <DST_DIR>/*.nc             (복사된 avg 파일들)

SRC_DIR="${1:?Need SRC_DIR}"
DST_DIR="${2:?Need DST_DIR}"
PREFIX="${3:-nifs5km}"

mkdir -p "$DST_DIR"

# avg 파일 존재 확인
shopt -s nullglob
avg_files=("$SRC_DIR/${PREFIX}_avg_"*.nc)
if [[ ${#avg_files[@]} -eq 0 ]]; then
  echo "No avg files found: $SRC_DIR/${PREFIX}_avg_*.nc" >&2
  exit 1
fi

# cycle(XXXX) 목록 추출/정렬
mapfile -t cycles < <(
  ls -1 "$SRC_DIR/${PREFIX}_avg_"*.nc \
  | sed -E "s/.*_avg_([0-9]{4})_.*/\1/" \
  | sort -n | uniq
)

: > "$DST_DIR/avg_use_list.txt"

first="${cycles[0]}"

# 첫 cycle: 0001~0008
for n in 0001 0002 0003 0004 0005 0006 0007 0008; do
  f="$SRC_DIR/${PREFIX}_avg_${first}_${n}.nc"
  [[ -f "$f" ]] && echo "$(basename "$f")" >> "$DST_DIR/avg_use_list.txt"
done

# 이후 cycle: 0005~0008
for c in "${cycles[@]:1}"; do
  for n in 0005 0006 0007 0008; do
    f="$SRC_DIR/${PREFIX}_avg_${c}_${n}.nc"
    [[ -f "$f" ]] && echo "$(basename "$f")" >> "$DST_DIR/avg_use_list.txt"
  done
done

# 복사
while read -r base; do
  src="$SRC_DIR/$base"
  dst="$DST_DIR/$base"
  cp -v "$src" "$dst"
done < "$DST_DIR/avg_use_list.txt"

echo "Wrote list: $DST_DIR/avg_use_list.txt"
