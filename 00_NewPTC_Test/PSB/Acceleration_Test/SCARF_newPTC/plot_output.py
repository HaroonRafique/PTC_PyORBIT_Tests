import matplotlib as mpl
mpl.use('Agg')
import os
import numpy as np
import scipy.io as sio
import pylab as plt
import glob

png_dir = 'png'
for d in [png_dir]:
	if not os.path.exists(d):
		os.mkdir(d)

plt.close('all')

filename = 'output/output.mat'
data = sio.loadmat(filename, squeeze_me=True)

gamma = data['gamma']
beta = np.sqrt(1.-1./gamma**2)
C = 50*np.pi
c = 3e8
time = np.cumsum(C/(beta*c)) * 1e6

f, axs = plt.subplots(2,figsize=(6,6), sharex=True)
ax = axs[0]
# ax.plot(data['turn'], 1e6*data['epsn_x']/gamma/beta)
# ax.plot(data['turn'], 1e6*data['epsn_y']/gamma/beta)
# ax.set_ylabel('physical emittance (um)')
ax.plot(data['turn'], 1e6*data['epsn_x'], 'b',  label='x')
ax.plot(data['turn'], 1e6*data['epsn_y'], 'g', label='y')
ax.plot(data['turn'], 1e6*(data['epsn_x']+data['epsn_y'])/2, 'k--', label='(x+y)/2')
# ax.plot(time, 1e6*data['epsn_x'], label='x')
# ax.plot(time, 1e6*data['epsn_y'], label='y')
# ax.set_xlabel('time (us)')
ax.set_ylabel('normalized emittance (um)')
ax.legend(loc=4, ncol=3, fontsize=11)
ax = axs[1]
ax.plot(data['turn'], data['intensity'], 'r')
total_intensity = np.copy(data['intensity'])
try:
	ax.plot(data['turn'], data['intensity_lost'], 'r--')
	total_intensity += data['intensity_lost']
except:
	pass
# try:
# 	ax.plot(data['turn'], data['foilhits']/100, 'b')
# except:
# 	pass
ax.set_xlabel('turn')
ax.set_ylabel('intensity')
ylim = np.round(10*max(total_intensity))/10
ax.set_ylim(min([ylim-0.025e13, data['intensity'][-1]-0.025e13]),ylim+0.025e13)
plt.tight_layout()
plt.savefig('%s/emittance_evolution.png'%png_dir, dpi=400)

try:
	f, ax = plt.subplots(1, figsize=(6,4))
	ax.plot(data['turn'], data['bunching_factor'])
	ax.set_ylim(0,0.8)
	ax.set_xlabel('turn')
	ax.set_ylabel('bunching factor')
	ax.grid()
	plt.tight_layout()
	plt.savefig('%s/BunchingFactor_evolution.png'%png_dir, dpi=400)

	beta = np.sqrt(1-1./data['gamma']**2)
	f, ax = plt.subplots(1, figsize=(6,4))
	ax.plot(1./(beta*data['gamma']**2 * data['bunching_factor']))
	ax.grid()
	ax.set_xlabel('turn')
	ax.set_ylabel('space charge factor')
	plt.tight_layout()
	plt.savefig('%s/spacecharge_factor.png'%png_dir)

except:
	pass
	
plt.show()


