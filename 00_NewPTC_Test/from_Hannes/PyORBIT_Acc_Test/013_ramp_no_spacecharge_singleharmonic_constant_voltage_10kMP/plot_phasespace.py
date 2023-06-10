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
	bins=250 #500
	my_cmap = plt.cm.jet
	my_cmap.set_under('w',0.1)
	f, axs = plt.subplots(2,2,figsize=(10,8))
	ax = axs[0,0]
	ax.hist2d(x*1000., xp*1000.,bins=bins, cmap=my_cmap, vmin=0.1)
	ax.set_xlabel('x [mm]')
	ax.set_ylabel('xp [mrad]')
	ax = axs[0,1]
	ax.hist2d(y*1000., yp*1000.,bins=bins, cmap=my_cmap, vmin=0.1)
	ax.set_xlabel('y [mm]')
	ax.set_ylabel('yp [mrad]')
	ax = axs[1,0]
	ax.hist2d(x*1000., y*1000.,bins=bins, cmap=my_cmap, vmin=0.1)
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
	plt.savefig('%s/phasespace_%s.png'%(png_dir, turn), dpi=400)
	plt.close(f)

plt.show()
print 'DONE'
