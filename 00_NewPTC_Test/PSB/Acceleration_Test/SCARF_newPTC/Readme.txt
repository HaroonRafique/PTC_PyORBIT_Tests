# clone optics repository for PSB by executing
. _clone_optics_repository.sh

# make modifications to PSB lattice in file 
# PSB/madx/psb_injection_for_pyOrbit.madx
# (e.g. adjust the tunes, assign errors, ...)

# adapt simulation parameters in 
# ./simulation_parameters.py

# prepare the simulation files by executing (see
# comments in the file)
. _prepare_simulation.sh

