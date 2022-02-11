import matplotlib as mpl
mpl.use('Agg')
import os
import numpy as np
import scipy.io as sio
import pylab as plt
import glob
from lib.TuneDiagram.tune_diagram import resonance_lines

png_dir = 'png'
for d in [png_dir]:
	if not os.path.exists(d):
		os.mkdir(d)

plt.close('all')

# source_dir = 'input/'
source_dir = 'output/'

filename = 'mainbunch'
files = glob.glob(source_dir + filename + '*.mat')
files.sort()

for i, file in enumerate(files[::-1]):
	print file
	try: 
		turn = int(file.split('mainbunch_')[-1][:-4])
		turn = '%05d'%turn
	except:
		turn = ''
	particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']
	x  = particles.x
	xp = particles.xp
	y  = particles.y
	yp = particles.yp
	z  = particles.z
	dE = particles.dE

	fontsize=15

	qx = particles.ParticlePhaseAttributes[2,:]
	qy = particles.ParticlePhaseAttributes[3,:]
	qx[np.where(qx>0.5)] -= 1
	qy[np.where((qy>0.7) & (qx<0.25))] -= 1

	resonances = resonance_lines((3.6, 4.6),(3.7, 4.7),(1,2,3),16)
	fontsize=17
	f, ax = plt.subplots(1,figsize=(6,6))
	my_cmap = plt.cm.jet
	my_cmap.set_under('w',1)
	r = resonances
	ax.hist2d(4+qx, 4+qy, bins=1000, cmap=my_cmap, vmin=1, range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]]) #, norm=mcolors.PowerNorm(gamma))
	resonances.plot_resonance(f)
	ax.xaxis.label.set_size(fontsize)
	ax.yaxis.label.set_size(fontsize)
	ax.tick_params(labelsize=fontsize)
	plt.title('turn %s'%turn, fontsize=fontsize)
	plt.tight_layout()
	f.savefig('%s/tune_diagram_%s.png'%(png_dir, turn), dpi=500)
	plt.close(f)

plt.show()
print 'DONE'
