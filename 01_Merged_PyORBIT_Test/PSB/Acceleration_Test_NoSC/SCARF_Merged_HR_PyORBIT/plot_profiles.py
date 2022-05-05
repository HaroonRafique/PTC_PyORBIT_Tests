import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
import scipy.io as sio
import numpy as np
import glob
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import matplotlib.cm as cm
import os
import glob
# import mystyle as ms
fsize=16
# ms.mystyle_arial(fontsz=fsize, dist_tick_lab=10)

def _Gauss(x,x0,a,sigma):
	return a*exp(-(x-x0)**2/(2*sigma**2))

def GaussianFit(x, y):
	mean = sum(x*y)/sum(y)
	sigma = np.sqrt(sum(y*(x-mean)**2)/sum(y))
	amplitude = max(y)
	popt,pcov = curve_fit(_Gauss,x,y,p0=[mean,amplitude,sigma])
	amplitude_norm = popt[1]*np.sqrt(2*np.pi)/(x[1]-x[0]) * popt[2] / np.float(sum(y))
	# print amplitude_norm
	return popt, amplitude_norm

def _DoubleGauss(x,mu,ampl1,sigma1,ampl2,sigma2):
	# return abs(ampl1)*np.exp(-(x-mu)**2/(2*sigma1**2)) + abs(ampl2)*np.exp(-(x-mu)**2/(2*sigma2**2))
	return abs(ampl1)*np.exp(-(x-mu)**2/(2*sigma1**2)) + abs(ampl2)*np.exp(-(x-mu)**2/(2*sigma2**2))

def DoubleGaussianFit(x,y):
	try:
		p,_ = GaussianFit(x,y)
		popt,pcov = curve_fit(_DoubleGauss,x,y,p0=[p[0],p[1]*0.95,p[2],p[1]*0.05,2*p[2]]) #, bounds=([-1.,-1e9, 0., -1e9, 0.], [1., 1e9, 1., 1e9, 1.]))
		d = x[1]-x[0]
		s = np.float(sum(y))
		if abs(popt[1])>abs(popt[3]):
			A1 = abs(popt[1])
			sig1 = abs(popt[2])
			A2 = abs(popt[3])
			sig2 = abs(popt[4])
		else:
			A1 = abs(popt[3])
			sig1 = abs(popt[4])			
			A2 = abs(popt[1])
			sig2 = abs(popt[2])
		A1_norm = A1*np.sqrt(2*np.pi)/d * sig1 / s
		A2_norm = A2*np.sqrt(2*np.pi)/d * sig2 / s
		# print A1_norm, A2_norm, A1_norm + A2_norm
		return {'mu': popt[0], 'A1': A1, 'sig1': sig1, 'A2': A2, 'sig2': sig2, 'A1_norm': A1_norm, 'A2_norm': A2_norm, 'pcov': pcov, 'p': popt}
	except:
		return {k: np.nan for k in ['mu', 'A1', 'sig1', 'A2', 'sig2', 'pcov', 'p', 'A1_norm', 'A2_norm']}

def calculate_second_moment(x, y, y_thresh_rel=0.02):
    # y -= min(y)
    # y -= np.median(y)
    y = y/max(y)
    i_filtered = np.where(y>y_thresh_rel)
    x_f = x[i_filtered]
    y_f = y[i_filtered]
    W = sum(y_f)
    mu = np.sum(x_f*y_f)/W
    sig = np.sqrt(np.sum(y_f*(x_f-mu)**2)/W)
    return sig, mu, x_f, y_f



png_dir = 'png'
for d in [png_dir]:
	if not os.path.exists(d):
		os.mkdir(d)

plt.close('all')

positions, amplitudes = {'x':[], 'y': [], 'z': []}, {'x':[], 'y': [], 'z': []}
turns = []

directory = 'output'

lastprofiles_pos = {}
lastprofiles_ampl = {}

files = sorted(glob.glob('%s/bunch_profile_*.mat'%directory))
files_to_plot = files[::]
f, axs = plt.subplots(1,3,figsize=(12,4))
f2, ax2 = plt.subplots(1,figsize=(6,4))
colors = cm.rainbow(np.linspace(0, 1, len(files_to_plot)))
for filename, c in zip(files_to_plot, colors):
	print filename
	with open(filename, 'r') as fid:
		a = sio.loadmat(fid, squeeze_me=True, struct_as_record=False)
		
		for i, u in enumerate(['x', 'y', 'z'][:]):
			pos = a[u][0]
			ampl = np.copy(a[u][1]) / (pos[1]-pos[0])
			if u in ['x', 'y']:
				pos *= 1e3
			axs[i].plot(pos, ampl, c=c)
			lastprofiles_pos[u] = pos
			lastprofiles_ampl[u] = ampl
		turn = int(filename.split('_')[-1][:-4])
		ax2.plot(pos,ampl+len(files_to_plot)/7.*turn,'k',lw=0.5)
		
axs[0].set_xlabel('x (mm)')
axs[1].set_xlabel('y (mm)')
axs[2].set_xlabel('z (m)')
axs[0].set_ylabel('density (a.u.)')
axs[1].set_ylabel('density (a.u.)')
axs[2].set_ylabel('density (a.u.)')
f.tight_layout()
f.savefig('%s/profiles.png'%png_dir, dpi=400)

ax2.set_xlabel('z (m)')
ax2.set_ylabel('counts')
f2.savefig('%s/mountainrange.png'%png_dir, dpi=400)

f3, axs3 = plt.subplots(1,3,figsize=(12,4))
axs3[0].plot(lastprofiles_pos['x'], lastprofiles_ampl['x'], lw=2)
popt, _ = GaussianFit(lastprofiles_pos['x'], lastprofiles_ampl['x'])
xlim = axs3[0].get_xlim()
x = np.linspace(xlim[0], xlim[1], 1000)
axs3[0].plot(x, _Gauss(x,popt[0],popt[1],popt[2]), 'r')
axs3[0].set_xlabel('x (mm)')
axs3[0].set_ylabel('counts')
axs3[1].plot(lastprofiles_pos['y'], lastprofiles_ampl['y'], lw=2)
popt, _ = GaussianFit(lastprofiles_pos['y'], lastprofiles_ampl['y'])
ylim = axs3[1].get_xlim()
y = np.linspace(ylim[0], ylim[1], 1000)
axs3[1].plot(y, _Gauss(y,popt[0],popt[1],popt[2]), 'r')
axs3[1].set_xlabel('y (mm)')
axs3[1].set_ylabel('counts')
axs3[2].plot(lastprofiles_pos['z'], lastprofiles_ampl['z'], lw=2)
axs3[2].set_ylabel('counts')
axs3[2].set_xlabel('z (m)')
plt.tight_layout()
f3.savefig('%s/profiles_last.png'%png_dir, dpi=400)
plt.show()