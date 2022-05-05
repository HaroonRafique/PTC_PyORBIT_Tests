# 12.10.2017: implemented possibility for double Gaussian transverse distribution 
# 14.10.2017: added sigma from FWHM of dp/p profile

import math
import sys
from itertools import chain
import numpy as np
import csv
import random
import orbit_mpi

from bunch import Bunch
from orbit.injection.joho import JohoLongitudinal
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist2D, GaussDist2D, KVDist2D
from orbit.utils.consts import mass_proton, speed_of_light, pi
from DoubleRF import DoubleRF

import scipy.io as sio
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def _Gauss(x,x0,a,sigma):
	return a*exp(-(x-x0)**2/(2*sigma**2))

def _GaussianFit(x, y):
	mean = sum(x*y)/sum(y)
	sigma = np.sqrt(sum(y*(x-mean)**2)/sum(y))
	amplitude = max(y)
	popt,pcov = curve_fit(_Gauss,x,y,p0=[mean,amplitude,sigma])
	amplitude_norm = popt[1]*np.sqrt(2*np.pi)/(x[1]-x[0]) * popt[2] / np.float(sum(y))
	return popt, amplitude_norm

def _Gaussian_sigma_from_FWHM(x,y):
	from scipy.interpolate import UnivariateSpline
	spline = UnivariateSpline(x, y-np.max(y)/2, s=0)
	r1, r2 = spline.roots() 
	return (r2-r1)/2.3548




class LongitudinalSquareDistribution():

	def __init__(self, z_min, z_max, dE_max, m):
		self.z_min = z_min
		self.z_max = z_max
		self.dE_max = dE_max
		self.m = m
		self.dist = lambda z,m: (1-np.clip(z,0,1)**2)**(m-0.5)

	def getCoordinates(self, n_mp=1):
		dist = self.dist
		z_min = self.z_min
		z_max = self.z_max
		dE_max = self.dE_max
		m = self.m

		z = np.linspace(-z_max,z_max,100)
		U_ = []
		V_ = []
		W_ = []		
		while len(U_)<n_mp:
			u = np.random.uniform(z_min,z_max,n_mp)
			v = np.random.uniform(-1,1,n_mp)
			w = np.random.uniform(0,1,n_mp)
			d = dist(np.abs(v), m)
			mask = np.where(w < d)[0]
			U_.extend(u[mask])
			V_.extend(v[mask]*dE_max)
			W_.extend(w[mask])
		z_rand = np.array(U_[:n_mp])
		dE_rand = np.array(V_[:n_mp])
		return z_rand, dE_rand


def generate_initial_distribution_noLongPainting(parameters, Lattice=None, output_file='ParticleDistribution.in', outputFormat='pyOrbit',
								  summary_file='ParticleDistribution_summary.txt', summary_mat_file=None):
	assert outputFormat in ['Orbit', 'pyOrbit']
	p = parameters
	beta = p['beta']
	gamma = p['gamma']
	if Lattice:
		p['alphax0'] = Lattice.alphax0
		p['betax0']  = Lattice.betax0
		p['alphay0'] = Lattice.alphay0
		p['betay0']  = Lattice.betay0
		p['etax0']   = Lattice.etax0
		p['etapx0']  = Lattice.etapx0
		p['etay0']   = Lattice.etay0
		p['etapy0']  = Lattice.etapy0
		p['x0']      = Lattice.orbitx0
		p['xp0']     = Lattice.orbitpx0
		p['y0']      = Lattice.orbity0
		p['yp0']     = Lattice.orbitpy0
		p['gamma_transition'] = Lattice.gammaT
		p['circumference']    = Lattice.getLength()

	# building the distributions

	R = p['circumference']/2/np.pi
	h_main = np.atleast_1d(p['harmonic_number'])[0]
	m = p['LongitudinalJohoParameter']
	dE_rms = p['LongitudinalDistribution_dE_rms']
	dE_max = 2*dE_rms*np.sqrt((m+1)/2)
	z_min = p['LongitudinalDistribution_z_min']
	z_max = p['LongitudinalDistribution_z_max']

	Longitudinal_distribution = LongitudinalSquareDistribution(z_min, z_max, dE_max, m)
	z, dE = Longitudinal_distribution.getCoordinates(p['n_macroparticles'])
	dE += p['LongitudinalDistribution_dE_mean']
	phi = - z * h_main / R
	dpp = dE / (p['energy'] * beta**2 * 1.e-9)

	# transverse coordinates
	x,xp,y,yp = [],[],[],[]
	for epsn_x, epsn_y, intensity in zip(np.atleast_1d(p['epsn_x']), np.atleast_1d(p['epsn_y']), np.atleast_1d(p['intensity'])):
		# twiss containers
		twissX = TwissContainer(alpha = p['alphax0'], beta = p['betax0'], emittance = epsn_x / gamma / beta)
		twissY = TwissContainer(alpha = p['alphay0'], beta = p['betay0'], emittance = epsn_y / gamma / beta)

		Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=p['TransverseCut'])
		n_macroparticles_tmp = int(p['n_macroparticles']*(intensity/np.sum(p['intensity'])))
		Transverse_coords = np.array(map(lambda i: Transverse_distribution.getCoordinates(), xrange(n_macroparticles_tmp)))
		x.extend(Transverse_coords[:,0].tolist())
		xp.extend(Transverse_coords[:,1].tolist())
		y.extend(Transverse_coords[:,2].tolist())
		yp.extend(Transverse_coords[:,3].tolist())
	# in case x has not yet a length of n_macroparticles
	while len(x)<p['n_macroparticles']:
		Transverse_coords = Transverse_distribution.getCoordinates()
		x.append(Transverse_coords[0])
		xp.append(Transverse_coords[1])
		y.append(Transverse_coords[2])
		yp.append(Transverse_coords[3])
	x = np.array(x) + p['x0']  + dpp * p['etax0']
	xp = np.array(xp) + p['xp0'] + dpp * p['etapx0']
	y = np.array(y) + p['y0']  + dpp * p['etay0']
	yp = np.array(yp) + p['yp0'] + dpp * p['etapy0']

	# only the main CPU is actually writing its distribution to a file ...
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	if orbit_mpi.MPI_Comm_rank(comm) == 0:
		with open(output_file,"w") as fid:
			csv_writer = csv.writer(fid, delimiter=' ')
			if outputFormat == 'Orbit':
				x  *= 1000.
				xp *= 1000.
				y  *= 1000.
				yp *= 1000.
				map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]]), range(p['n_macroparticles']))	
			elif outputFormat == 'pyOrbit':
				map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]]), range(p['n_macroparticles']))	

		if summary_file:
			with open(summary_file, 'w') as fid:
				map(lambda key: fid.write(key + ' = ' + str(p[key]) + '\n'), p)

		if summary_mat_file:
			with open(summary_mat_file, 'w') as fid:
				sio.savemat(fid, parameters) 

		print '\nCreated particle distribution with ' + str(p['n_macroparticles']) + ' macroparticles into file: ', output_file
	
	orbit_mpi.MPI_Barrier(comm)

	return output_file



def generate_initial_distribution_LongPainting(parameters, Lattice=None, output_file='ParticleDistribution.in', outputFormat='pyOrbit',
								  summary_file='ParticleDistribution_summary.txt', summary_mat_file=None):
	assert outputFormat in ['Orbit', 'pyOrbit']
	p = parameters
	beta = p['beta']
	gamma = p['gamma']
	if Lattice:
		p['alphax0'] = Lattice.alphax0
		p['betax0']  = Lattice.betax0
		p['alphay0'] = Lattice.alphay0
		p['betay0']  = Lattice.betay0
		p['etax0']   = Lattice.etax0
		p['etapx0']  = Lattice.etapx0
		p['etay0']   = Lattice.etay0
		p['etapy0']  = Lattice.etapy0
		p['x0']      = Lattice.orbitx0
		p['xp0']     = Lattice.orbitpx0
		p['y0']      = Lattice.orbity0
		p['yp0']     = Lattice.orbitpy0
		p['gamma_transition'] = Lattice.gammaT
		p['circumference']    = Lattice.getLength()

	# building the distributions
	from scipy.constants import c, proton_mass
	from scipy.constants import physical_constants

	R = p['circumference']/2/np.pi
	h_main = np.atleast_1d(p['harmonic_number'])[0]
	m = p['LongitudinalJohoParameter']
	t2z = beta * c 
	Ekin = p['energy']
	E0 = physical_constants['proton mass energy equivalent in MeV'][0] * 1e6
	d = np.genfromtxt(p['painting_definition_file'], delimiter=',', names=True)
	duration_average = np.mean(d['dtMax'] - d['dtMin'])
	
	n_macroparticles_average = float(p['n_macroparticles'])/float(p['nturns_accumulation'])
	n_macroparticles_sum = 0
	for k, (dtMin, dtMax, dEMin, dEMax, EHalfHeight) in enumerate(zip(d['dtMin'], d['dtMax'], d['dEMin']/1e9, d['dEMax']/1e9, d['EHalfHeight']/1e9)):
		if k>=p['nturns_accumulation']:
			print "WARNING: number of accumulation turns smaller than turns provided by painting file ... additional turns ignored"
			continue
		output_file_tmp = ''.join(output_file.split('.')[:-1]) + '%i.%s'%(k+1, output_file.split('.')[-1])
		z_min = -dtMin*t2z
		z_max = -dtMax*t2z

		if k<p['nturns_accumulation']-1:
			n_macroparticles = int(np.round(n_macroparticles_average * (dtMax - dtMin) / duration_average))
		else: 
			n_macroparticles = p['n_macroparticles'] - n_macroparticles_sum
		n_macroparticles_sum += n_macroparticles	
		Longitudinal_distribution = LongitudinalSquareDistribution(z_min, z_max, EHalfHeight, m)
		z, dE = Longitudinal_distribution.getCoordinates(n_macroparticles)
		phi = - z * h_main / R
		dE += np.interp(z, [z_min, z_max], [dEMin, dEMax])
		dpp = dE / (p['energy'] * beta**2 * 1.e-9)

		# transverse coordinates
		x,xp,y,yp = [],[],[],[]
		for epsn_x, epsn_y, intensity in zip(np.atleast_1d(p['epsn_x']), np.atleast_1d(p['epsn_y']), np.atleast_1d(p['intensity'])):
			# twiss containers
			twissX = TwissContainer(alpha = p['alphax0'], beta = p['betax0'], emittance = epsn_x / gamma / beta)
			twissY = TwissContainer(alpha = p['alphay0'], beta = p['betay0'], emittance = epsn_y / gamma / beta)

			Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=p['TransverseCut'])
			n_macroparticles_tmp = int(n_macroparticles*(intensity/np.sum(p['intensity'])))
			Transverse_coords = np.array(map(lambda i: Transverse_distribution.getCoordinates(), xrange(n_macroparticles_tmp)))
			x.extend(Transverse_coords[:,0].tolist())
			xp.extend(Transverse_coords[:,1].tolist())
			y.extend(Transverse_coords[:,2].tolist())
			yp.extend(Transverse_coords[:,3].tolist())
		# in case x has not yet a length of n_macroparticles
		while len(x)<n_macroparticles:
			Transverse_coords = Transverse_distribution.getCoordinates()
			x.append(Transverse_coords[0])
			xp.append(Transverse_coords[1])
			y.append(Transverse_coords[2])
			yp.append(Transverse_coords[3])
		x = np.array(x) + p['x0']  + dpp * p['etax0']
		xp = np.array(xp) + p['xp0'] + dpp * p['etapx0']
		y = np.array(y) + p['y0']  + dpp * p['etay0']
		yp = np.array(yp) + p['yp0'] + dpp * p['etapy0']

		# only the main CPU is actually writing its distribution to a file ...
		comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
		if orbit_mpi.MPI_Comm_rank(comm) == 0:
			with open(output_file_tmp,"w") as fid:
				csv_writer = csv.writer(fid, delimiter=' ')
				if outputFormat == 'Orbit':
					x  *= 1000.
					xp *= 1000.
					y  *= 1000.
					yp *= 1000.
					map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]]), range(n_macroparticles))	
				elif outputFormat == 'pyOrbit':
					map(lambda i: csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]]), range(n_macroparticles))	

			if summary_file:
				with open(summary_file, 'w') as fid:
					map(lambda key: fid.write(key + ' = ' + str(p[key]) + '\n'), p)

			if summary_mat_file:
				with open(summary_mat_file, 'w') as fid:
					sio.savemat(fid, parameters) 

			print '\nCreated particle distribution with ' + str(n_macroparticles) + ' macroparticles into file: ', output_file_tmp
		
	print "\ntotal number of macroparticles generated: %i"%n_macroparticles_sum
	if k+1<p['nturns_accumulation']:
		print "WARNING: number of accumulation turns larger than turns provided by painting file !!!"
	print "average chopping factor: %1.3f"%(duration_average/(p['circumference']/beta/c))
	orbit_mpi.MPI_Barrier(comm)

	return output_file
