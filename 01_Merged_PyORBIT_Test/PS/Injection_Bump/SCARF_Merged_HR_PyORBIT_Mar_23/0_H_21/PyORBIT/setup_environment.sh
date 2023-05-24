# script to setup the PyOrbit environment 
# execute like: . setup_environment.sh

pyOrbit_dir=/apps20/contrib/ptc-pyorbit/PTC-PyORBIT/Merged_HR_PyORBIT_2/py-orbit

source ${pyOrbit_dir}/customEnvironment.sh
echo "customEnvironment done"
source ${pyOrbit_dir}/../virtualenvs/py2.7/bin/activate
echo "python packages charged"
#source ${pyOrbit_dir}/../setup_ifort.sh
#echo "ifort charged (necessary for running)"

which python

ORBIT_ROOT_fullpath=`readlink -f ${ORBIT_ROOT}` 
echo 
echo "*****************************************************"
echo 
echo "full PyOrbit path:  ${ORBIT_ROOT_fullpath}"
echo
. ${ORBIT_ROOT}/../CheckGitStatus.sh ${ORBIT_ROOT_fullpath}

# Test: remove creation of bytecode to avoid race condition bug
export PYTHONDONTWRITEBYTECODE=1
