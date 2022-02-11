#!/bin/bash
#SBATCH --job-name=OP_Acc2
#SBATCH --output=slurm.%N.%j.out
#SBATCH --error=slurm.%N.%j.err
#SBATCH --constraint="intel"
#SBATCH --ntasks=80
#SBATCH --partition=ibis
#SBATCH --time=168:00:00
#SBATCH --exclusive

ml ifort/2018.3.222-GCC-7.3.0-2.30
BATCH_ROOT_DIR=/home/vol05/scarf1015
RUN_DIR=/work3/isis/scarf1015/PTC_PyORBIT_Tests/00_NewPTC_Test/PSB/Acceleration_Test/SCARF_newPTC
OrigIwd=$(pwd)

# Make an output folder in the root directory to hold SLURM info file
cd ${BATCH_ROOT_DIR}
output_dir="Simulation_SLURM_Outputs"
mkdir -p $output_dir

# Fill the SLURM info file
simulation_info_file="${BATCH_ROOT_DIR}/${output_dir}/simulation_info_${SLURM_JOB_ID}.${SLURM_NODEID}.${SLURM_PROCID}.txt"
echo "PyOrbit path:  `readlink -f ${ORBIT_ROOT}`" >> ${simulation_info_file}
echo "Run path:  `readlink -f ${RUN_DIR}`" >> ${simulation_info_file}
echo "Submit host:  `readlink -f ${SLURM_SUBMIT_HOST}`" >> ${simulation_info_file}
echo "SLURM Job name:  `readlink -f ${SLURM_JOB_NAME}`" >> ${simulation_info_file}
echo "SLURM Job ID:  `readlink -f ${SLURM_JOB_ID}`" >> ${simulation_info_file}
echo "SLURM Nodes allocated:  `readlink -f ${SLURM_JOB_NUM_NODES}`" >> ${simulation_info_file}
echo "SLURM CPUS per Node:  `readlink -f ${SLURM_CPUS_ON_NODE}`" >> ${simulation_info_file}
echo "SLURM Node ID:  `readlink -f ${SLURM_NODEID}`" >> ${simulation_info_file}
echo "SLURM total cores for job:  `readlink -f ${SLURM_NTASKS}`" >> ${simulation_info_file}
echo "SLURM process ID:  `readlink -f ${SLURM_PROCID}`" >> ${simulation_info_file}
echo "****************************************" >> ${simulation_info_file}

# Enter job directory, clean it, and setup environment -> SLURM info file
cd ${RUN_DIR}
./clean_all.sh
. setup_environment.sh >> ${simulation_info_file}

tstart=$(date +%s)

# Run the job
mpirun ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/pyOrbit.py

tend=$(date +%s)
dt=$(($tend - $tstart))
echo "total simulation time (s): " $dt >> ${simulation_info_file}