#!/usr/bin/env python
# Python script to create a SLURM submission script for PyORBIT
# 21 March 2019 Haroon Rafique CERN BE-ABP-HSI 
# 19 March 2021 Haroon Rafique STFC ISIS

import os
  
#-----------------------------------------------------------------------
#	SETTINGS
#-----------------------------------------------------------------------  
script_name = "SLURM_submission_script.sh"

# Switches
hyperthreading = True  # Enable hyperthreading
exclusive = True        # Exclusive nodes (see SLURM documentation)
autotime = True         # 7 days True, 2 days False
autotask = False        # Automatically set nodes to maximum tasks
clean_all = True        # Clean simulation folder before running (False when resuming pickle checkpoint)
constrain = True        # Force SCARF hardware
mem_limits = False      # Specify memory limits, default 4GB

# Must be chosen

# ~ queue = 'scarf' # Default 16/20/24 cores/Node
queue = 'ibis' # isis exclusive 10x SCARF16 nodes 20 cores/Node

# ~ n_nodes = 4 
n_tasks_tot = 80

# ~ space_charge_flag = int(os.getcwd().split('/')[-1][0])
# ~ print 'simulation_parameters: space charge = ', space_charge_flag
# ~ transverse_plane = os.getcwd().split('/')[-1][2]
# ~ print 'simulation_parameters: transverse_plane = ', transverse_plane
# ~ scan_tune = os.getcwd().split('/')[-1][-2:]

# ~ jobname = str(space_charge_flag) + '_' + str(transverse_plane) + '_' + str(scan_tune)
jobname = 'PS_Acc'

path_to_simulation = os.path.dirname(os.path.realpath(__file__)) # This directory

# Optional - have to use with correct switches
manual_time = '168:00:00'       # manually set using format 'hours:minutes:seconds'
manual_tasks = 80               # manually change ntasks
manual_mem_limit = '4000'       # Manually set memory limit
# ~ manual_constraints = '"[scarf17|scarf18]"'     # select multiple SCARF clusters
# ~ manual_constraints = 'scarf17'                 # select single SCARF cluster
manual_constraints = '"intel"'                 # select intel hardware

# Defaults - can be changed
output_file_name = 'slurm.%N.%j.out'
error_file_name = 'slurm.%N.%j.err'
root_dir = '/home/vol05/scarf1015'
simulation_file = 'pyOrbit.py'
#-----------------------------------------------------------------------
#	AUTOMATICALLY FORMAT SCRIPT
#-----------------------------------------------------------------------  
n_tasks = 0
if autotask:
        if 'scarf' in queue: 
                n_tasks = 16 # 20, 24 ??
        elif 'ibis' in queue: 
                n_tasks = 80
        else: 
                print('queue not recognised')
                exit(0)
else: n_tasks = manual_tasks

time = '168:00:00'
if autotime:
        time = '168:00:00'
else: time = manual_time

constraints = ''
if constrain:
        constraints = manual_constraints
else:
        if queue == 'scarf': constraints = ''
        elif queue == 'ibis': constraints = 'scarf16'

mem_limit = '4000'
if mem_limits:
        if 'scarf15' or 'scarf16' in constraints: mem_limit = '6000'
        elif 'scarf17' in constraints: mem_limit = '5000'
        elif 'scarf18' or 'scarf19' in constraints: mem_limit = '7994' 
        elif 'gpu13' in constraints: mem_limit = '4000' 
        elif 'scarf14' or 'gpu14' or 'gpu15' or 'gpu17' or 'scarf20' in constraints: mem_limit = '8000' 
else: mem_limit = manual_mem_limit

#-----------------------------------------------------------------------
#	WRITE FILE
#-----------------------------------------------------------------------  
if os.path.exists(script_name):  
	print ('SLURM submission script ' + script_name + ' already exists. Deleting')
	os.remove(script_name)

print ("Creating ", script_name)

f= open(script_name,"w")

f.write('#!/bin/bash')
f.write('\n#SBATCH --job-name=' + str(jobname))
f.write('\n#SBATCH --output=' + str(output_file_name))
f.write('\n#SBATCH --error=' + str(error_file_name))
if constrain: f.write('\n#SBATCH --constraint=' + str(constraints))
f.write('\n#SBATCH --ntasks=' + str(n_tasks_tot))
# ~ f.write('\n#SBATCH --nodes=' + str(n_nodes))
# ~ f.write('\n#SBATCH --ntasks-per-node=' + str(n_tasks))
f.write('\n#SBATCH --partition=' + str(queue))
f.write('\n#SBATCH --time=' + str(time))
if mem_limits:  f.write('\n#SBATCH --mem-per-cpu='+str(mem_limit)+'M')
if exclusive: f.write('\n#SBATCH --exclusive')
#if not hyperthreading: f.write('\n#SBATCH --hint=nomultithread')
if not hyperthreading: f.write('\nexport OMP_NUM_THREADS=1')
f.write('\n')
f.write('\nml ifort/2018.3.222-GCC-7.3.0-2.30')
f.write('\nBATCH_ROOT_DIR=' + str(root_dir))
f.write('\nRUN_DIR=' + str(path_to_simulation))
f.write('\nOrigIwd=$(pwd)')
f.write('\n')
f.write('\n# Make an output folder in the root directory to hold SLURM info file')
f.write('\ncd ${BATCH_ROOT_DIR}')
f.write('\noutput_dir="Simulation_SLURM_Outputs"')
f.write('\nmkdir -p $output_dir')
f.write('\n')
f.write('\n# Fill the SLURM info file')
f.write('\nsimulation_info_file="${BATCH_ROOT_DIR}/${output_dir}/simulation_info_${SLURM_JOB_ID}.${SLURM_NODEID}.${SLURM_PROCID}.txt"')
f.write('\necho "PyOrbit path:  `readlink -f ${ORBIT_ROOT}`" >> ${simulation_info_file}')
f.write('\necho "Run path:  `readlink -f ${RUN_DIR}`" >> ${simulation_info_file}')
f.write('\necho "Submit host:  `readlink -f ${SLURM_SUBMIT_HOST}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Job name:  `readlink -f ${SLURM_JOB_NAME}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Job ID:  `readlink -f ${SLURM_JOB_ID}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Nodes allocated:  `readlink -f ${SLURM_JOB_NUM_NODES}`" >> ${simulation_info_file}')
f.write('\necho "SLURM CPUS per Node:  `readlink -f ${SLURM_CPUS_ON_NODE}`" >> ${simulation_info_file}')
f.write('\necho "SLURM Node ID:  `readlink -f ${SLURM_NODEID}`" >> ${simulation_info_file}')
f.write('\necho "SLURM total cores for job:  `readlink -f ${SLURM_NTASKS}`" >> ${simulation_info_file}')
f.write('\necho "SLURM process ID:  `readlink -f ${SLURM_PROCID}`" >> ${simulation_info_file}')
f.write('\necho "****************************************" >> ${simulation_info_file}')
f.write('\n')
f.write('\n# Enter job directory, clean it, and setup environment -> SLURM info file')
f.write('\ncd ${RUN_DIR}')
if clean_all:f.write('\n./clean_all.sh')
f.write('\n. setup_environment.sh >> ${simulation_info_file}')
# ~ f.write('\n')
# ~ f.write('\n# Load correct MPI')
#f.write('\nmodule load mpi/mpich-3.0-x86_64')
# ~ f.write('\nmodule load mpi/mpich-x86_64')
# ~ f.write('\nmodule load mpi/openmpi-x86_64')
f.write('\n')
f.write('\ntstart=$(date +%s)')
f.write('\n')
f.write('\n# Run the job')
if hyperthreading:f.write('\nmpirun ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/' + str(simulation_file))	
# ~ if hyperthreading:f.write('\nsrun ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/' + str(simulation_file))	
else:f.write('\nmpirun ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/' + str(simulation_file))
# ~ else:f.write('\nsrun --hint=nomultithread ${ORBIT_ROOT}/bin/pyORBIT ${RUN_DIR}/' + str(simulation_file))
f.write('\n')
f.write('\ntend=$(date +%s)')
f.write('\ndt=$(($tend - $tstart))')
f.write('\necho "total simulation time (s): " $dt >> ${simulation_info_file}')

f.close()

print ('SLURM submission script creation finished')
