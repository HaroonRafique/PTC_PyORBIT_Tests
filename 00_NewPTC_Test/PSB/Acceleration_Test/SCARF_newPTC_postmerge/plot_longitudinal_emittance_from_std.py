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

# f, ax = plt.subplots(1)
# ax.plot(data['bunchlength_rms'])
# ax.plot(data['dE_rms'])
# ax.plot(data['bunchlength_rms']/beta*data['dE_rms'])

f, ax = plt.subplots(1)
# ax.plot(data['turn'], data['eps_z'])
ax.plot(4*np.pi*data['bunchlength_rms']/(beta*3e8)*data['dE_rms']*1e9)
ax.set_xlabel('turn')
ax.set_ylabel('lognitudinal emittance (eVs)')
plt.savefig('%s/longitudinal_emittance_from_std.png'%png_dir, dpi=400)

plt.show()