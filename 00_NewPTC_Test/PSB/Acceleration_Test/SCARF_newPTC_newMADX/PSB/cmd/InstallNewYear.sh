#!/bin/bash
# Bash shell script to make a new version of the PSB optics from the layoutDB
#
# (c) I.Efthymiopoulos - March 2018
#
# Usage :
#     1/ download thew new file from the layoutDB
#     2/ run the python script SplitLaoutDB.py to create the different input files
#     3/ run this script to install the new files in the sequence, elements and strength directories. Use the option -keep to save a copy of the existing files before copying
#
#   ./InstallNewYeear.sh -f psb_2018  -keep
#

usage () { echo "Usage : $0 [-f FILE] [-k 1/0 to keep a copy of existing files]" 1>&2; exit 0; }

while getopts f:k: option ; do
  case "${option}" in
    f) FILENAM=${OPTARG};;
    k) KEEP=${OPTARG};;
    h | *)
      usage
      exit 0;
      ;;

  esac
done
shift $((OPTIND-1))

if [ -z "${FILENAM}" ] || [ -z "${KEEP}" ]; then
  usage
fi

echo "Input parameters:"
echo "File name : ${FILENAM}"
echo "Keep existing files = ${KEEP}"

if [ ${KEEP}==1 ]; then
  now=$(date +"%Y.%m.%d.%I.%M.%p")
  cp ../sequence/psb.seq ../sequence/psb_${now}.seq
  cp ../elements/psb.ele ../elements/psb_${now}.ele
  cp ../strength/psb_orbit.str ../strength/psb_orbit_${now}.str
fi


cp ./${FILENAM}.seq ../sequence/psb0.seq

# -- replace LAYOUT_ID with SLOT_ID in the sequence files
sed 's/LAYOUT_ID/SLOTID/g' ../sequence/psb0.seq > ../sequence/psb.seq
rm ../sequence/psb0.seq

cp ./${FILENAM}.ele ../elements/psb.ele
cp ./${FILENAM}.str ../strength/psb_orbit.str
