
#!/usr/bin/env bash
# 사용: bash rerun_failed.sh run_all.log [rerun.log]
set -euo pipefail

SRC_LOG="${1:-run_all.log}"
RERUN_LOG="${2:-rerun.log}"
[[ -f "$SRC_LOG" ]] || { echo "로그 파일 없음: $SRC_LOG" >&2; exit 2; }
: > "$RERUN_LOG"

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

# (rc=...) 가 있는 [FAILD] 라인만 추출 → 요약부의 [FAILD]는 제외
mapfile -t CMDS < <( awk '
  /^\[FAILD\] / && /\(rc=/ {
    line=$0
    sub(/^\[FAILD\] /,"",line)       # 프리픽스 제거
    sub(/ \(rc=.*$/,"",line)         # 뒤쪽 (rc=...) 제거
    if(!seen[line]++){ print line }  # 중복 제거
  }
' "$SRC_LOG" )

if ((${#CMDS[@]}==0)); then
  echo "재실행할 실패 커맨드가 없습니다: $SRC_LOG"
  exit 0
fi

FAILED=()
for cmd in "${CMDS[@]}"; do
  echo "[$(ts)] RERUN $cmd" | tee -a "$RERUN_LOG"
  set +e
  bash -lc "$cmd" 2>&1 | tee -a "$RERUN_LOG"
  rc=${PIPESTATUS[0]}
  set -e
  if [[ $rc -eq 0 ]]; then
    echo "[OK] $cmd" | tee -a "$RERUN_LOG"
  else
    echo "[FAILD] $cmd (rc=$rc)" | tee -a "$RERUN_LOG"
    FAILED+=("$cmd")
  fi
done

echo "----- RERUN SUMMARY -----" | tee -a "$RERUN_LOG"
if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "ALL [OK]" | tee -a "$RERUN_LOG"
  exit 0
else
  for c in "${FAILED[@]}"; do echo "[FAILD] $c" | tee -a "$RERUN_LOG"; done
  exit 1
fi
