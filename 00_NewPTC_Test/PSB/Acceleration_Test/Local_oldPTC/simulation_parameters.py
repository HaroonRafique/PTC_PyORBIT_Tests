import numpy as np

choppingFactor = 1.0 #0.62 #0.7
m = 1.5 # parameter of binomial longitudinal distribution 
epsn_x = 1e-6 
epsn_y = 1e-6 
TransverseCut = 15
intensity = 2.0e12
nturns_accumulation = 1
n_macroparticles = 10000
orbitx0_injection = 0
orbity0_injection = 0
macrosize = intensity/float(n_macroparticles)
# ~ turns_max = 250000
turns_max = 250
# ~ turns_print = range(-1, turns_max, 1000)
turns_print = range(-1, turns_max, 1)
turns_injection = np.arange(nturns_accumulation)

particledistribution_directory = 'Distribution_at_injection'

parameters = {
	'LongitudinalJohoParameter': m,
	'n_macroparticles': n_macroparticles,
	'intensity': intensity,
	'nturns_accumulation': nturns_accumulation,
	'turns_injection': turns_injection,
	'orbitx0_injection': orbitx0_injection,
	'orbity0_injection': orbity0_injection,
	'epsn_x': epsn_x,
	'epsn_y': epsn_y,
	'TransverseCut': TransverseCut,
	'macrosize': macrosize,
	'turns_max': turns_max,
	'turns_print': turns_print,
	'particledistribution_directory': particledistribution_directory,
}


# use this for simulations without longitudinal painting
########################################################
z_lim = choppingFactor*2*np.pi*25
z_min, z_max,  = -0.5*25*np.pi, 0.5*25*np.pi
dE_rms = 0.1e-3 # max from Linac4 is 0.5e-3 GeV
dE_mean = 0

parameters.update({
	'LongitudinalDistribution_z_max': z_max,
	'LongitudinalDistribution_z_min': z_min,
	'LongitudinalDistribution_dE_rms': dE_rms,
	'LongitudinalDistribution_dE_mean': dE_mean,
})

