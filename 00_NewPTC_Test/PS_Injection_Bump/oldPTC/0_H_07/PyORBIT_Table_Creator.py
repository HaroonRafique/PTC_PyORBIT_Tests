# Script to read single or double columned MADX tfs tables containing 
# magnet ramp data, and converting this to a PTC-PyORBIT readable PTC 
# table. Note that we ignore the timing column as in MAD-X this is not 
# defined as a time but as a step. We therefore create our own timing.
#
# Haroon Rafique CERN BE-ABP-HSI 28.06.19
#
import numpy as np
import os
from math import log10, floor
from scipy.misc import factorial

# From Hannes Bartosik CERN BE-ABP-HSI
def write_PTCtable(filename, multipole_orders, time, normal_components, skew_components):
	print '\n Create table ', filename
	multipole_orders = np.atleast_1d(multipole_orders)
	factors = 1./factorial(multipole_orders-1) # the factorial factor is needed to be consistent with MADX
	if factors.shape[0] <= 1:
		normal_components = (factors * (normal_components))
		skew_components   = (factors * (skew_components))
		arr = np.empty((normal_components.shape[0], 3), dtype=normal_components.dtype)
		arr[:,0] = time
		arr[:,1] = normal_components
		arr[:,2] = skew_components
		n_lines = len(time)
		n_multipoles = len(multipole_orders) # number of multipole orders to be changed

	else:
		normal_components = (factors.T * np.atleast_2d(normal_components))
		skew_components   = (factors.T * np.atleast_2d(skew_components))
		arr = np.empty((normal_components.shape[0], 1+normal_components.shape[1]*2), dtype=normal_components.dtype)
		arr[:,0] = time
		arr[:,1::2] = normal_components
		arr[:,2::2] = skew_components
		n_lines = len(time)
		n_multipoles = len(multipole_orders) # number of multipole orders to be changed
		
	with open(filename, 'w') as fid:
		fid.write('%d  1  %d\n'%(n_lines, n_multipoles))
		fid.write(' '.join(map(lambda i: '%d'%i, multipole_orders)) + '\n')
		for j in xrange(n_lines):
			fid.write('\t'.join(map(lambda i: '%+1.11f'%i, arr[j, :]))+'\n')
	return

# Read TFS table to obtain our dipole kick or quadrupole gradient
# We ignore the time column
def Read_Single_Column_TFS_Return_Data(file_in):
	print '\n Reading single column file ', file_in
	fi = open(file_in, 'r')
	contents = fi.readlines()	
	x=[]
	data =np.ndarray(shape=(1,1))
	data = np.delete(data, 0, 0)

	for l in contents:
		if ('@' in l) or ('$' in l) or ('*' in l):
			pass
		else:
			x = [float(l.split()[1])]
			data = np.vstack([data, x])

	data = np.delete(data, 0, 0)
	return data

# Read TFS table to obtain our dipole kick or quadrupole gradient
# We ignore the time column
def Read_Double_Column_TFS_Return_Data(file_in, f1=1., f2=1.):
	print '\n Reading double column file', file_in
	fi = open(file_in, 'r')
	contents = fi.readlines()
	x=[]
	data =np.ndarray(shape=(2,2))
	data = np.delete(data, 1, 0)
	
	for l in contents:
		if ('@' in l) or ('$' in l) or ('*' in l):
			pass
		else:
			x = [f1*float(l.split()[1]), f2*float(l.split()[2])]
			data = np.vstack([data, x])

	data = np.delete(data, 0, 0)
	return data

# Function takes the data, creates the correct timing for each value in
# seconds, then appends the required number of intervals with zeroes
def Create_Timing(ramp_stop_time, simulation_stop_time, data):
	print '\n Create timing'
	data_2 = data

	# create sequence for time
	d_len = len(data)
	try:
		d_shape = data_2.shape[1]
	except AttributeError as error:
		d_shape = 1
	
	time = np.linspace(0, ramp_stop_time, d_len)

	# Calculate interval to complete data until the end of the simulation
	interval = ramp_stop_time / (d_len-1)

	# Append zeroes (bump is closed)
	steps = int((simulation_stop_time - ramp_stop_time)	/ interval) + 5

	for i in xrange(steps):
		if i is 0:
			pass
		else:
			time = np.append(time, (ramp_stop_time + (i*interval)))
			if d_shape == 1:
				x = [0.0]
				data_2 = np.vstack([data_2, x])
			elif d_shape == 2:
				x = [0.0, 0.0]
				data_2 = np.vstack([data_2, x])

	result = np.column_stack([time, data_2])

	return result

# Create the input files for the closure (2nd half) of the injection bump

# Table time is in seconds
# 1 turn = 2.287E-6 seconds
# We want a full bump of 2E-3 seconds, or a half bump of 1E-3 seconds
# 2200 turns = 5.0314E-3 s
# Bump on for 1E-3 seconds
# zeroes from 1E-3 - 5.0314E-3 s
# 500 turns = 1.1435E-3 s

l_bsw = 0.19

B_40 = Read_Double_Column_TFS_Return_Data('./MADX_Tables/BSEXT40.tfs', f1=1./l_bsw, f2=1./l_bsw)
B_40_final = Create_Timing(1.1435E-3, 5.0314E-3, B_40)
# write_PTCtable(file, multipoles, time column, normal components, skew components
write_PTCtable('./PTC-PyORBIT_Tables/BSEXT40.dat', (1,3), B_40_final[:,0],  B_40_final[:,[1,2]], B_40_final[:,[1,2]]*0)

B_42 = Read_Double_Column_TFS_Return_Data('./MADX_Tables/BSEXT42.tfs', f1=1./l_bsw, f2=1./l_bsw)
B_42_final = Create_Timing(1.1435E-3, 5.0314E-3, B_42)
write_PTCtable('./PTC-PyORBIT_Tables/BSEXT42.dat', (1,3), B_42_final[:,0], B_42_final[:,[1,2]], B_42_final[:,[1,2]]*0)

B_43 = Read_Double_Column_TFS_Return_Data('./MADX_Tables/BSEXT43.tfs', f1=1./l_bsw, f2=1./l_bsw)
B_43_final = Create_Timing(1.1435E-3, 5.0314E-3, B_43)
write_PTCtable('./PTC-PyORBIT_Tables/BSEXT43.dat', (1,3), B_43_final[:,0], B_43_final[:,[1,2]], B_43_final[:,[1,2]]*0)

B_44 = Read_Double_Column_TFS_Return_Data('./MADX_Tables/BSEXT44.tfs', f1=1./l_bsw, f2=1./l_bsw)
B_44_final = Create_Timing(1.1435E-3, 5.0314E-3, B_44)
write_PTCtable('./PTC-PyORBIT_Tables/BSEXT44.dat', (1,3), B_44_final[:,0], B_44_final[:,[1,2]], B_44_final[:,[1,2]]*0)
