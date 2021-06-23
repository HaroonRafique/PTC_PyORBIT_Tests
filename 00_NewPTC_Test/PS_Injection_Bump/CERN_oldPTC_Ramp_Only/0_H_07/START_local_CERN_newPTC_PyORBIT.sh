#!/bin/bash

set -o errexit

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the SC script> <N CPUs>"
    exit $E_BADARGS
fi

if [ ! -n "$2" ]
  then
    echo "Usage: `basename $0` <name of the SC script> <N CPUs>"
    exit $E_BADARGS
fi

source /home/HR/Codes/PTC_PyORBIT/NewPTC_Test2/venv/bin/activate
echo "python packages charged"
source /home/HR/Codes/PTC_PyORBIT/NewPTC_Test2/py-orbit_newPTC/customEnvironment.sh
echo "customEnvironment done"

# Check if link to $ORBIT_ROOT/lib exists
my_link=./lib
if [ -L ${my_link} ] ; then
   if [ -e ${my_link} ] ; then
      echo "lib $my_link: exists and linked to $ORBIT_ROOT"
   else
      echo "lib $my_link: broken link,if using local lib/ please be aware \
   of possible linking errors to user defined python scripts. If this \
   happens move your user defined python lib files to $ORBIT_ROOT and \
   link using:$> ln -s $ORBIT_ROOT/lib/ lib"
   fi
elif [ -e ${my_link} ] ; then
   echo "lib $my_link: Not a link, if using local lib/ please be aware \
   of possible linking errors to user defined python scripts. If this \
   happens move your user defined python lib files to $ORBIT_ROOT and \
   link using:$> ln -s $ORBIT_ROOT/lib/ lib"
else
   echo "lib $my_link: Missing, creating link to $ORBIT_ROOT"
   ln -s $ORBIT_ROOT/lib/ lib
fi

/usr/lib64/openmpi/bin/mpirun -np $2 ${ORBIT_ROOT}/bin/pyORBIT $1
