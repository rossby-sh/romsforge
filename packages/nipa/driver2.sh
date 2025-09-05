#!/usr/bin/env bash
# get_* 은 FETCH env, mk_* 는 PROC env로 실행
# [OK]/[FAILD] 로깅, 실패해도 계속 진행, 마지막에 요약

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT" || exit 2

# ====== 여기만 네 환경에 맞게 설정 ======
# 방법 A) 각 env의 python 절대경로 지정
#EXEC_FETCH="/path/to/fetch-env/bin/python"
#EXEC_PROC="/path/to/proc-env/bin/python"

# 방법 B) conda run 사용 (A 대신 사용 가능)
EXEC_FETCH="conda run -n fetch python"
EXEC_PROC="conda run -n romsforge python"

# 공통/개별 옵션 필요하면 여기에
ARGS_FETCH=""   # 예: "-c config_all.yaml"
ARGS_PROC=""    # 예: "-c config_all.yaml"

LOG="${1:-run_all.log}"
: > "$LOG"
TS(){ date +"%Y-%m-%d %H:%M:%S"; }

SCRIPTS=(
  "get_cds_era5.py"
  "get_hycom.py"
  "get_cmems_bio.py"
  "get_cmems_ostia.py"
  "get_modis_chl.py"
  "mk_ini_nifs.py"
  "mk_bry_nifs.py"
  "mk_frc_nifs.py"
)

FAILED=()

for s in "${SCRIPTS[@]}"; do
  if [[ "$s" == get_* ]]; then
    runner="$EXEC_FETCH"
    extra="$ARGS_FETCH"
  else
    runner="$EXEC_PROC"
    extra="$ARGS_PROC"
  fi

  cmd="$runner $s $extra"
  echo "[$(TS)] EXEC $cmd" | tee -a "$LOG"

  set +e
  bash -lc "$cmd" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  set -e

  if [[ $rc -eq 0 ]]; then
    echo "[OK] $cmd" | tee -a "$LOG"
  else
    echo "[FAILD] $cmd (rc=$rc)" | tee -a "$LOG"
    FAILED+=("$cmd")
  fi
done

echo "----- SUMMARY -----" | tee -a "$LOG"
if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "ALL [OK]" | tee -a "$LOG"
  exit 0
else
  for f in "${FAILED[@]}"; do echo "[FAILD] $f" | tee -a "$LOG"; done
  exit 1
fi
