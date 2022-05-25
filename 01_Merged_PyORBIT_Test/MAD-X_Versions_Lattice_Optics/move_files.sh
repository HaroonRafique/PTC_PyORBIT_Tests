#!/bin/bash

# Check if directories needed exist, if not, create them

if [ ! -d "MADX_Tables" ]; then
  mkdir MADX_Tables
fi
if [ ! -d "MADX_Twiss" ]; then
  mkdir MADX_Twiss
fi
if [ ! -d "PTC_Twiss" ]; then
  mkdir PTC_Twiss
fi
if [ ! -d "PTC-PyORBIT_Tables" ]; then
  mkdir PTC-PyORBIT_Tables
fi

mv BSEXT40.tfs MADX_Tables
mv BSEXT42.tfs MADX_Tables
mv BSEXT43.tfs MADX_Tables
mv BSEXT44.tfs MADX_Tables
mv BSEXT_Out.tfs MADX_Tables
mv *.tfs MADX_Twiss
mv MADX_Twiss/madx_twiss.tfs .
mv *.ptc PTC_Twiss
rm PTC-PyORBIT_flat_file.flt
mv MADX_Twiss/optimised_flat_file.tfs .
mv MADX_Twiss/optimised_bare_simplified.tfs .
