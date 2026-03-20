#!/usr/bin/env bash
set -euo pipefail

# ==========
# paths
# ==========
CONDA_HOME="/usr/local/miniconda3"

BASE="/home/shjo/github/romsforge/packages/nipa_auto"
readonly PREP="${BASE}/src/preprocess"

NL_MDL="${BASE}/roms_application/roms_820_nemuscs_latest"

readonly INPUT_ROOT="${BASE}/roms_inputs"
readonly OUTPUT_ROOT="${BASE}/roms_outputs"

readonly SRC_TPL="${PREP}/config.template"
readonly SRC_OUT="${PREP}/config.yaml"

readonly NL_MDL_TPL="${NL_MDL}/nl_nipa.template"
readonly NL_MDL_OUT="${NL_MDL}/nl_nipa.in"

readonly ROMS_EXEC=${NL_MDL}/execM_nl
# ==========
# prevent overlap
# ==========
mkdir -p "${BASE}"
exec 9>"${BASE}/.lock_main"
flock -n 9 || { echo "[INFO] already running. exit."; exit 0; }

# ---- check running ROMS job ----
if qstat -u "$USER" | grep -q "nipa_"; then
    echo "[INFO] ROMS job already running. skip."
    exit 0
fi

# ==========
# run_date
# ==========
run_date="${run_date:-$(date)}"
run_date="2026-03-12"

export TZ=UTC

readonly RUN_DATE="$(date -u -d "$run_date" +%Y%m%d)"
readonly RUN_DATE_ISO="$(date -u -d "$run_date" +%F)"
TIME_REF="days since 2000-01-01 00:00:00"

# ==========
# directories
# ==========
readonly INPUTS_DIR="${INPUT_ROOT}/${RUN_DATE}"
readonly OUTPUTS_DIR="${OUTPUT_ROOT}/${RUN_DATE}"
readonly LOG_DIR="${INPUTS_DIR}"

mkdir -p "${INPUTS_DIR}" "${OUTPUTS_DIR}" #"${LOG_DIR}"

readonly LOGFILE="${LOG_DIR}/system_log_${RUN_DATE}.md"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "[INFO] RUN_DATE=${RUN_DATE}"
echo "[INFO] INPUTS_DIR=${INPUTS_DIR}"
echo "[INFO] OUTPUTS_DIR=${OUTPUTS_DIR}"

# ==========
# forcing window
# ==========
t0="$(date -u -d "${RUN_DATE_ISO} 00:00:00" +%s)"

readonly INI_DATE="$(date -u -d "@$t0" '+%F %T')"
readonly RUN_DATE_START="$(date -u -d "@$((t0 - 3*3600))" '+%F %T')"
readonly RUN_DATE_END="$(date -u -d "@$((t0 + (24*2+3)*3600))" '+%F %T')"
echo "[INFO] INI_DATE=${INI_DATE}"


# ==========
# Define infile variables
# ==========
TITLE_NL="NIPA test run (RUN_DATE: ${RUN_DATE})"
VARINFO_NL=${NL_MDL}/varinfo_nipa.dat

NTILEI=4
NTILEJ=4

DT_NL=200
readonly REF_TIME="${TIME_REF#days since }"

ref_sec=$(date -u -d "${REF_TIME}" +%s)
run_sec=$(date -u -d "${INI_DATE}" +%s)
readonly TIME_REF_MDL=$(date -u -d "$REF_TIME" +%Y%m%d)
readonly DSTART=$(echo "scale=0; (${run_sec} - ${ref_sec}) / 86400" | bc)

readonly NRST_NL=$(( (86400 / DT_NL) * 1 ))
readonly NRST_OL=$NRST_NL

save_interval=4

readonly NHIS_NL=$(( (86400 / DT_NL) / save_interval ))
readonly NDEFHIS_NL=${NRST_NL}
readonly NAVG_NL=${NRST_NL}
readonly NDEFAVG_NL=${NRST_NL}
readonly NINFO_NL=${NHIS_NL}

GRDNAME_NL=${INPUT_ROOT}/fixed/NWP4_grd_314_10m.nc
readonly ININAME_NL=${INPUTS_DIR}/roms_nipa_ini_${RUN_DATE}.nc
readonly BRYNAME_NL=${INPUTS_DIR}/roms_nipa_bry_${RUN_DATE}.nc
readonly FRCNAME_NL_01=${INPUTS_DIR}/roms_nipa_frc_${RUN_DATE}.nc
readonly FRCNAME_NL_02=${INPUT_ROOT}/fixed/Frc_lon_lat_36525.nc
readonly SSFNAME_NL=${INPUT_ROOT}/fixed/nifs01_24km_river.nc
readonly TIDENAME_NL=${INPUT_ROOT}/fixed/nifs01_24km_tide_20000101.nc

readonly RSTNAME_NL=${OUTPUTS_DIR}/nipa_rst_${RUN_DATE}.nc
readonly HISNAME_NL=${OUTPUTS_DIR}/nipa_his_${RUN_DATE}.nc
readonly AVGNAME_NL=${OUTPUTS_DIR}/nipa_avg_${RUN_DATE}.nc

#=====================================================================
#if USE_RST == TRUE: --> 모델 돌아가는 거 보고 RST를 따로 저장할지 생각 
#                        RST관련 규칙을 만들어야 할 듯
#    NRRECT=-1  
#    RAN_DATE=${RUNDATE}-1 
#    readonly ININAME_NL=${INPUTS_DIR}/roms_nipa_ini_${RAN_DATE}.nc
#else:
#    readonly ININAME_NL=${INPUTS_DIR}/roms_nipa_ini_${RUN_DATE}.nc
#=====================================================================

BASE_DIR=${BASE}/data

WGHTNAME_NL=${INPUT_ROOT}/fixed/hycom_wght_24km_314.nc
WGHTNAME_BIO_NL=${INPUT_ROOT}/fixed/cmems_wght_24km_314.nc
CALC_WGHT=True

OGCM_NAME=${BASE}/data/hycom/hycom_korea_${RUN_DATE}.nc
OGCM_BIO_NAME=${BASE}/data/cmems/cmems_bio_${RUN_DATE}.nc

BIO_TYPE=npzd

# ---  export variables & convert configuration templates ---

ENV_CONF_NL_VARS="
GRDNAME_NL ININAME_NL BRYNAME_NL 
GRDNAME_NL ININAME_NL BRYNAME_NL CALC_WGHT 
FRCNAME_NL_01 WGHTNAME_NL WGHTNAME_BIO_NL 
OGCM_NAME OGCM_BIO_NAME BIO_TYPE BASE_DIR 
INI_DATE TIME_REF RUN_DATE_START RUN_DATE_END
"

export $ENV_CONF_NL_VARS

ENV_CONF_NL_TPL=$(printf '$%s ' $ENV_CONF_NL_VARS)

envsubst "${ENV_CONF_NL_TPL}" < "${SRC_TPL}" > "${SRC_OUT}"

# --- export variables & convert model templates ---

ENV_NL_VARS="
TITLE_NL VARINFO_NL NTILEI NTILEJ
REF_TIME DSTART TIME_REF_MDL DT_NL
NRST_NL NRST_OL NHIS_NL NDEFHIS_NL
NAVG_NL NDEFAVG_NL NINFO_NL
GRDNAME_NL ININAME_NL BRYNAME_NL
FRCNAME_NL_01 FRCNAME_NL_02
SSFNAME_NL TIDENAME_NL
RSTNAME_NL HISNAME_NL AVGNAME_NL
"

export $ENV_NL_VARS

ENV_NL_TPL=$(printf '$%s ' $ENV_NL_VARS)

envsubst "$ENV_NL_TPL" < "$NL_MDL_TPL" > "$NL_MDL_OUT"


#============================================================================












