#!/usr/bin/env bash
set -euo pipefail

# ==========
# paths
# ==========
CONDA_HOME="/usr/local/miniconda3"

BASE="/home/shjo/github/romsforge/packages/nipa_auto"
SRC="${BASE}/src"
PREP="${SRC}/preprocess"

MDL="${BASE}/roms_application/roms_820_nemuscs_latest"

INPUT_ROOT="${BASE}/roms_inputs"
OUTPUT_ROOT="${BASE}/roms_outputs"

SRC_TPL="${PREP}/config.template"
SRC_OUT="${PREP}/config.yaml"

MDL_TPL="${MDL}/nl_nipa.template"
MDL_OUT="${MDL}/nl_nipa.in"

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
# date
# ==========
export TZ=UTC

RUN_DATE="$(date -u +%Y%m%d)"
RUN_DATE_ISO="$(date -u +%F)"

TIME_REF="days since 2000-01-01 00:00:00"

# ==========
# directories
# ==========
INPUTS_DIR="${INPUT_ROOT}/${RUN_DATE}"
OUTPUTS_DIR="${OUTPUT_ROOT}/${RUN_DATE}"
LOG_DIR="${INPUTS_DIR}"

mkdir -p "${INPUTS_DIR}" "${OUTPUTS_DIR}" #"${LOG_DIR}"

LOGFILE="${LOG_DIR}/system_log_${RUN_DATE}.md"
exec > >(tee -a "${LOGFILE}") 2>&1

echo "[INFO] RUN_DATE=${RUN_DATE}"
echo "[INFO] INPUTS_DIR=${INPUTS_DIR}"
echo "[INFO] OUTPUTS_DIR=${OUTPUTS_DIR}"

# ==========
# forcing window
# ==========
t0="$(date -u -d "${RUN_DATE_ISO} 00:00:00" +%s)"

INI_DATE="$(date -u -d "@$t0" '+%F %T')"
RUN_DATE_START="$(date -u -d "@$((t0 - 3*3600))" '+%F %T')"
RUN_DATE_END="$(date -u -d "@$((t0 + 27*3600))" '+%F %T')"

echo "[INFO] INI_DATE=${INI_DATE}"

# ==========
# ROMS dstart
# ==========
REF_TIME="2000-01-01 00:00:00"

ref_sec=$(date -u -d "${REF_TIME}" +%s)
run_sec=$(date -u -d "${INI_DATE}" +%s)

DAY=$(echo "scale=1; (${run_sec} - ${ref_sec}) / 86400" | bc)

# ==========
# pipeline status
# ==========
PIPE_FILE="${INPUTS_DIR}/pipeline_status.md"

init_pipeline() {

if [ -f "$PIPE_FILE" ]; then
    echo "[INFO] pipeline_status exists → resume"
    return
fi

cat <<EOF > "$PIPE_FILE"
[ ] Download HYCOM frcst
[ ] CMEMS download
[ ] INI build
[ ] INI bio add
[ ] BRY build
[ ] BRY bio add
EOF

}

is_done() {
grep -q "^\[X\] $1" "$PIPE_FILE"
}

mark_done() {
sed -i "s/\[ \] $1/[X] $1/" "$PIPE_FILE"
sed -i "s/\[FAILED\] $1/[X] $1/" "$PIPE_FILE"
}

mark_fail() {
sed -i "s/\[ \] $1/[FAILED] $1/" "$PIPE_FILE"
}

run_stage() {

local name="$1"
shift

if is_done "$name"; then
    echo "[SKIP] $name already done"
    return
fi

echo "[INFO] --- START: $name ---"

if "$@"; then
    echo "[OK] --- DONE: $name ---"
    mark_done "$name"
else
    echo "[ERR] --- FAILED: $name ---"
    mark_fail "$name"
fi

}

# ==========
# init pipeline
# ==========
init_pipeline

# ==========
# config.template → config.yaml
# ==========
export BASE INPUT_ROOT INPUTS_DIR OUTPUTS_DIR RUN_DATE TIME_REF RUN_DATE_START RUN_DATE_END INI_DATE DAY

ENV_VARS='
$RUN_DATE_START
$RUN_DATE_END
$BASE
$INPUT_ROOT
$OUTPUTS_DIR
$INPUTS_DIR
$RUN_DATE
$TIME_REF
$INI_DATE
$DAY
'

envsubst "${ENV_VARS}" < "${SRC_TPL}" > "${SRC_OUT}"
envsubst "${ENV_VARS}" < "${MDL_TPL}" > "${MDL_OUT}"

# ==========
# conda runners
# ==========
FETCH_RUN="$CONDA_HOME/bin/conda run --no-capture-output -n fetch python"
ROMS_RUN="$CONDA_HOME/bin/conda run --no-capture-output -n romsforge python"


# ==================================================
# PREPROCESS
# ==================================================

# ==========
# forcing
# ==========

run_stage "Download HYCOM frcst" \
  $FETCH_RUN "${PREP}/get_hycom_frcst.py"

run_stage "CMEMS download" \
$FETCH_RUN "${PREP}/get_cmems_bio_nipa.py"

# ==========
# preprocess
# ==========
run_stage "INI build" \
$ROMS_RUN "${PREP}/mk_ini_single_new.py"

run_stage "INI bio add" \
$ROMS_RUN "${PREP}/add_bio_ini3.py"

run_stage "BRY build" \
$ROMS_RUN "${PREP}/mk_bry_single_new2.py"

run_stage "BRY bio add" \
$ROMS_RUN "${PREP}/add_bio_bry.py"

# ==================================================
# Preprocessing fallback
# ==================================================











# ==================================================
# MODEL RUN
# ==================================================


# ==========
# PBS script
# ==========
PBS_SCRIPT="${OUTPUTS_DIR}/run_roms_${RUN_DATE}.pbs"

cat <<EOF > "${PBS_SCRIPT}"
#!/bin/bash
set -eu

#PBS -q workq
#PBS -l nodes=nft00:ppn=36,walltime=100:00:00
#PBS -N nipa_${RUN_DATE}
#PBS -j oe
#PBS -o ${OUTPUTS_DIR}/pbs_${RUN_DATE}.log

module load intel21/compiler-21
module load intel21/openmpi-5.0.0
module load intel21/netcdf-4.6.1

export HDF5_USE_FILE_LOCKING=FALSE
export HDF5_COLL_METADATA_WRITE=0

NP=\$(wc -l < "\$PBS_NODEFILE")

cd ${MDL}

mpirun -np \$NP \\
  --map-by core \\
  --bind-to core \\
  --report-bindings \\
  --mca coll ^hcoll,han \\
  ./execM_nl nl_nipa.in \\
  > ${OUTPUTS_DIR}/mdl_log_${RUN_DATE}.md

EOF

# ==========
# submit
# ==========
JOBID=$(qsub "${PBS_SCRIPT}")

echo "[INFO] PBS submitted"
echo "[INFO] JOBID=${JOBID}"
