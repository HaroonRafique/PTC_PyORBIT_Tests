#!/bin/bash

# preparation for space charge simulations in the PSB with 
# the new lattice configuration post LS2 as extracted from the 
# layout database. the following commands have to be executed 
# before launching an actual pyOrbit simulation. The main 
# simulation settings are defined in ./simulation_parameters.py

python 001_run_MADX.py
# runs MADX on the file PSB/madx/psb_injection_for_pyOrbit.madx to 
# match the tunes (as defined in the file itself) and the dynamic 
# betabeat correction during the fall of the injection chicane. 
# Errors have to be assigned to the sequence just before generating
# the flat file.

./START_local_CERN_PTC_PyORBIT.sh 002_pyOrbit_prepare_RFtables.py 1
# generate RF tables and dnergy ramp based on Simon's input files

./START_local_CERN_PTC_PyORBIT.sh 003_pyOrbit_generateDistribution.py 1
# generates input distrubition without longitudinal painting, 
# for longitudinal painting see below


