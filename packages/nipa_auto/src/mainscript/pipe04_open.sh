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
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -N nipa_${RUN_DATE}
#PBS -j oe
#PBS -o ${OUTPUTS_DIR}/pbs_${RUN_DATE}.log

module load intel21/compiler-21
module load intel21/openmpi-5.0.0
module load intel21/netcdf-4.6.1

export HDF5_USE_FILE_LOCKING=FALSE
export HDF5_COLL_METADATA_WRITE=0

NP=\$(wc -l < "\$PBS_NODEFILE")
cd ${NL_MDL}

mpirun -np 16 \
  -x OMP_NUM_THREADS \
  -x HDF5_USE_FILE_LOCKING \
  -x HDF5_COLL_METADATA_WRITE \
  --map-by slot \
  --bind-to core \
  --report-bindings \
  --mca coll ^hcoll,han \
  ${ROMS_EXEC} ${NL_MDL_OUT} > ${OUTPUTS_DIR}/mdl_log_${RUN_DATE}.md 2>&1

EOF

# ==========
# submit
# ==========
JOBID=$(qsub "${PBS_SCRIPT}")

echo "[INFO] PBS submitted"
echo "[INFO] JOBID=${JOBID}"
