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

filename = 'lostbunch'
files = glob.glob(source_dir + filename + '*.mat')
files.sort()

circumference = 157.08 # 50*np.pi

y_lim = 1

for i, file in enumerate(files[::-1]):
	print file
	try: 
		turn = int(file.split('bunch_')[-1][:-4])
		turn = '%05d'%turn
	except:
		turn = ''
	particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']
	x  = np.atleast_1d(particles.x) / 1e10
	xp = np.atleast_1d(particles.xp) / 1e10
	y  = np.atleast_1d(particles.y) / 1e10
	yp = np.atleast_1d(particles.yp) / 1e10
	z  = np.atleast_1d(particles.z)
	dE = np.atleast_1d(particles.dE)

	if not len(x):
		print '... no particles found ...'
		continue

	fontsize=15
	bins=100
	my_cmap = plt.cm.jet
	my_cmap.set_under('w',0.1)
	f, axs = plt.subplots(2,2,figsize=(10,8))
	ax = axs[0,0]
	ax.hist2d(x*1e3, xp*1e3,bins=bins, cmap=my_cmap, vmin=0.1)
	ax.set_xlabel('x [mm]')
	ax.set_ylabel('xp [mrad]')
	ax = axs[0,1]
	ax.hist2d(y*1e3, yp*1e3,bins=bins, cmap=my_cmap, vmin=0.1)
	ax.set_xlabel('y [mm]')
	ax.set_ylabel('yp [mrad]')
	ax = axs[1,0]
	ax.hist2d(x*1e3, y*1e3,bins=bins, cmap=my_cmap, vmin=0.1)
	ax.set_xlabel('x [mm]')
	ax.set_ylabel('y [mm] ')
	ax = axs[1,1]
	ax.hist2d(z, dE,bins=bins, cmap=my_cmap, vmin=0.1)
	ax.set_xlabel('z [m]')
	ax.set_ylabel('dE [GeV] ')
	for ax in axs.flatten():
		ax.xaxis.label.set_size(fontsize)
		ax.yaxis.label.set_size(fontsize)
		ax.tick_params(labelsize=fontsize)
	plt.suptitle('turn %s'%turn, fontsize=fontsize)
	plt.tight_layout(rect=[0, 0.03, 1, 0.95])
	plt.savefig('%s/lostparticles_phasespace_%s.png'%(png_dir, turn), dpi=400)

	f, ax = plt.subplots(1,figsize=(8,5))
	plt.hist(np.atleast_1d(particles.LostParticleAttributes)%(circumference),500, range=(0,circumference));
	ax.set_xlabel('s [m]')
	ax.set_ylabel('counts')
	ax.xaxis.label.set_size(fontsize)
	ax.yaxis.label.set_size(fontsize)
	y_lim = max([y_lim, ax.get_ylim()[1]])
	ax.set_ylim(0,y_lim)
	ax.tick_params(labelsize=fontsize)
	ax.set_xlim(0,circumference)
	plt.tight_layout()
	plt.savefig('%s/lostparticlelocations_%s.png'%(png_dir, turn), dpi=400)

	turn = int(turn)
	if not i:
		f, ax = plt.subplots(1,figsize=(8,5))
		ax.hist(np.atleast_1d(particles.LostParticleAttributes)/circumference+turn,(turn+1)/10,color='r', linewidth=0)
		ax.set_xlabel('turn')
		ax.set_ylabel('macro particles loss')
		ax.xaxis.label.set_size(fontsize)
		ax.yaxis.label.set_size(fontsize)
		ax.tick_params(labelsize=fontsize)
		ax.set_xlim(left=0)
		plt.tight_layout()
		plt.savefig('%s/loss_evolution.png'%png_dir, dpi=400)
	plt.close('all')

plt.show()
print 'DONE'
