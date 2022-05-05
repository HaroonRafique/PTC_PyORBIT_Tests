#!/bin/bash

# Pickle and output files
rm -r output
rm -r bunch_output
rm -r lost
rm -r input

rm *.tfs
rm *.ptc
rm *.png
#~ rm PTC-PyORBIT_flat_file.flt
rm tunespread.dat
rm madx.ps
rm All_Twiss/*

. clean_junk.sh
