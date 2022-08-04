import os
import numpy as np

# If we are using the folder naming convention:
# Space charge flag - horizontal/vertical scan - scan point
# e.g. 0_H_07

space_charge_flag = int(os.getcwd().split('/')[-2][0])
transverse_plane = os.getcwd().split('/')[-2][2]
scan_tune = os.getcwd().split('/')[-2][-2:]
print( 'simulation_parameters: space charge = ', space_charge_flag)
print( 'simulation_parameters: transverse_plane = ', transverse_plane)
print( 'simulation_parameters: scan_tune = ', scan_tune)
