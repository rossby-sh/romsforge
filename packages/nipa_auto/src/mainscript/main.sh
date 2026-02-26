
#!/usr/bin/env bash
set -euo pipefail

# ==========
# paths
# ==========
BASE="/home/shjo/github/romsforge/packages/nipa_auto"
SRC="${BASE}/src"
MAIN="${SRC}/mainscript"
PREP="${SRC}/preprocess"

INPUT_ROOT="${BASE}/roms_inputs"
OUTPUT_ROOT="${BASE}/roms_outputs"

TPL="${PREP}/config.template"
OUT="${PREP}/config.yaml"

# ==========
# prevent overlap (cron 중복 실행 방지)
# ==========
mkdir -p "${BASE}"
exec 9>"${BASE}/.lock_main"
flock -n 9 || { echo "[INFO] already running. exit."; exit 0; }

# ==========
# date
# ==========
RUN_DATE="$(date +%Y%m%d)"
RUN_DATE_ISO="$(date +%F)"

# ==========
# TIME_REF (ROMS reference time)
# - 추천: 파일로 관리 (없으면 여기 상수로 바꿔)
# ==========
TIME_REF="days since 2000-1-1" 

# ==========
# per-day dirs
# ==========
INPUTS_DIR="${INPUT_ROOT}/${RUN_DATE}"
OUTPUT_DIR="${OUTPUT_ROOT}/${RUN_DATE}"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${INPUTS_DIR}" "${OUTPUT_DIR}" "${LOG_DIR}"

LOGFILE="${LOG_DIR}/system_log_${RUN_DATE}.log"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "[INFO] RUN_DATE=${RUN_DATE}"
echo "[INFO] INPUTS_DIR=${INPUTS_DIR}"
echo "[INFO] OUTPUT_DIR=${OUTPUT_DIR}"
echo "[INFO] TIME_REF=${TIME_REF}"

# ==========
# forcing/bry 기간: 하루 구간 + 전후 3시간 패딩
# ==========
export TZ=UTC

RUN_DATE_ISO="$(date +%F)"
t0="$(date -d "${RUN_DATE_ISO} 00:00:00" +%s)"


INI_DATE="$(date -d "@$t0" '+%F %T')"

RUN_DATE_START="$(date -d "@$((t0 - 3*3600))" '+%F %T')"
RUN_DATE_END="$(date -d "@$((t0 + 27*3600))" '+%F %T')"
echo "[DBG] TZ=$TZ"
echo "[DBG] date=$(date '+%F %T %Z')"
echo "[DBG] date -u=$(date -u '+%F %T %Z')"
echo "[INFO] RUN_DATE_START=${RUN_DATE_START}"
echo "[INFO] RUN_DATE_END=${RUN_DATE_END}"

# ==========
# config.template -> config.yaml (화이트리스트 치환)
# ==========
export BASE INPUT_ROOT INPUTS_DIR RUN_DATE TIME_REF RUN_DATE_START RUN_DATE_END INI_DATE

ENV_VARS='
$RUN_DATE_START
$RUN_DATE_END
$BASE
$INPUT_ROOT
$INPUTS_DIR
$RUN_DATE
$TIME_REF
$INI_DATE
'

envsubst "${ENV_VARS}" < "${TPL}" > "${OUT}"

# 미치환 변수 검사
if grep -Eq '\$\{[^}]+\}' "${OUT}"; then
  echo "[ERR] unresolved variables in ${OUT}"
  grep -nE '\$\{[^}]+\}' "${OUT}" || true
  exit 1
fi

echo "[OK] generated ${OUT}"

# ==========
# run stages (너 설계 순서)
# ==========
#bash "${MAIN}/download_data.sh"
#bash "${MAIN}/preprocess.sh"
#bash "${MAIN}/run_roms.sh"
#bash "${MAIN}/postprocess.sh"

echo "[DONE] ${RUN_DATE}"
