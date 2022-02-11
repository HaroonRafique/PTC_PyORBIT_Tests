# script to setup the PyOrbit environment 
# execute like: . setup_environment.sh
pyOrbit_dir=/apps/contrib/ptc-pyorbit/PTC-PyORBIT/CERN_PTC-PyORBIT/py-orbit

# HR 18.03.21
# Note that the custom environment is python 2.7
# python 3.0 commands will not work inside a PyORBIT script but can be
# used for post processing

source ${pyOrbit_dir}/../virtualenvs/py2.7/bin/activate
echo "python packages charged"

source ${pyOrbit_dir}/customEnvironment.sh
echo "customEnvironment done"

which python

# The intel fortran licences are required for PTC. If you are running
# PyORBIT locally you will likely have to VPN into the STFC network
# to access the license server
#source ${pyOrbit_dir}/../setup_ifort.sh
#ml ifort/2018.3.222-GCC-7.3.0-2.30
#echo "ifort charged (necessary for running)"

ORBIT_ROOT_fullpath=`readlink -f ${ORBIT_ROOT}` 
echo 
echo "*****************************************************"
echo 
echo "full PyOrbit path:  ${ORBIT_ROOT_fullpath}"
echo
. ${ORBIT_ROOT}/../CheckGitStatus.sh ${ORBIT_ROOT_fullpath}

