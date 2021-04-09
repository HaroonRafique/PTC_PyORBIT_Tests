#!/bin/bash
. clean_run.sh
. clean_junk.sh
cd ../
. clean_folder.sh
rm *.png
rm *.tfs
cd PyORBIT
