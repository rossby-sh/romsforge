#!/bin/bash
set -eu

#PBS -q workq
#PBS -l nodes=nft00:ppn=16,walltime=100:00:00
#PBS -N nipa_20260312
#PBS -j oe
#PBS -o /home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312/pbs_20260312.log

module purge
module load intel21/compiler-21
module load intel21/openmpi-5.0.0
module load intel21/netcdf-4.6.1

export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export HDF5_COLL_METADATA_WRITE=0

ROMS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_application/roms_820_nemuscs_latest
ROMS_EXEC=${ROMS_DIR}/execM_nl
ROMS_IN=${ROMS_DIR}/nl_nipa.in
MPIRUN=/usr/local/mpi/intel21/openmpi-5.0.0/bin/mpirun

cd "${ROMS_DIR}"

echo "===== DEBUG INFO ====="
echo "HOSTNAME          : $(hostname)"
echo "PBS_JOBID         : ${PBS_JOBID:-undefined}"
echo "PBS_NODEFILE      : ${PBS_NODEFILE:-undefined}"
echo "PBS_TASKNUM       : ${PBS_TASKNUM:-undefined}"
echo "which mpirun      : $(which mpirun)"
echo "nodefile content  :"
cat "${PBS_NODEFILE}" || true
echo "======================"

"${MPIRUN}" -np 16 \
  --host "$(hostname)" \
  --map-by slot:OVERSUBSCRIBE \
  --bind-to core \
  --report-bindings \
  -x OMP_NUM_THREADS \
  -x HDF5_USE_FILE_LOCKING \
  -x HDF5_COLL_METADATA_WRITE \
  --mca coll ^hcoll,han \
  "${ROMS_EXEC}" "${ROMS_IN}" \
  > /home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312/mdl_log_20260312.md 2>&1
