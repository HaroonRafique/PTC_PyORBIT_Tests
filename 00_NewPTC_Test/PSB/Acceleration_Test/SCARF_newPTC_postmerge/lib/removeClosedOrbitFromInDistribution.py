import numpy as np
import scipy.io as sio
import glob 

remove_orbit_X = -81
remove_orbit_Y = 10

infolder = 'Distribution_at_injection_full/'
outfolder = 'Distribution_at_injection_full_no_BSW_noLP/'

files = glob.glob(infolder + 'Orbit*')

for fn in files:
	data =  np.genfromtxt(fn)
	data[:,0] -= remove_orbit_X
	data[:,2] -= remove_orbit_Y

	np.savetxt(outfolder+fn.split('/')[-1],data) 
