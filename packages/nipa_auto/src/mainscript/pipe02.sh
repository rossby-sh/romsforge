

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

