import os
import numpy as np

parameters = {}

parameters['turns_max'] 		= int(3)
parameters['n_macroparticles']	= int(1E5)
parameters['InitialDistnSigma'] = 10.	# Poincare distribution limits in units of sigma

parameters['intensity']			= 1E+13
parameters['bunch_length']		= 140e-9
parameters['epsn_x']			= 1E-3# 400 * np.pi * 1E-6 / 6. # 99% emittance 400 pi mm mrad -> Assume to be 6 sigma wide
parameters['epsn_y']			= 1E-3#400 * np.pi * 1E-6 / 6.
parameters['dpp_rms']			= 1e-03
parameters['LongitudinalJohoParameter'] = 1.2
parameters['LongitudinalCut'] 	        = 5.
parameters['TransverseCut']		= 5.
parameters['rf_voltage']		= 0.0
parameters['circumference']		= 163.748686344000
parameters['phi_s']			    = 0
parameters['macrosize']			= parameters['intensity']/float(parameters['n_macroparticles'])

# ISIS-II Injection 0.4 GeV
parameters['gamma']				= 1.426322879896
parameters['beta'] 	= np.sqrt(parameters['gamma']**2-1)/parameters['gamma']
print 'beta = ', parameters['beta'] 
c 			= 299792458
parameters['sig_z'] 	= (parameters['beta'] * c * parameters['bunch_length'])/4.

# Define how often we dump bunch output files
#-----------------------------------------------------------------------
n_turns = 1
parameters['turns_print'] = range(-1, parameters['turns_max'], n_turns) # every n_turns
parameters['turns_update'] = range(-1, parameters['turns_max'], n_turns) # every n_turns

# Simulation switches
#-----------------------------------------------------------------------
switches = {
    'Gaussian': True, # Gaussian or Poincare (False) distribution
	'Horizontal': True,   # Initial Poincare distribution in horizontal (true) / vertical (false)
	'Update_Twiss':	False   # Perform PTC twiss and dump each turn - needed to output tune changes
}

# PTC RF Table Parameters
#-----------------------------------------------------------------------
harmonic_factors = [1] # this times the base harmonic defines the RF harmonics (for SPS = 4620, PS 10MHz 7, 8, or 9)
time = np.array([0,1,2])
ones = np.ones_like(time)
Ekin_GeV = 0.4*ones
RF_voltage_MV = np.array([parameters['rf_voltage']*ones]).T # in MV
RF_phase = np.array([np.pi*ones]).T

RFparameters = {
	'harmonic_factors': harmonic_factors,
	'time': time,
	'Ekin_GeV': Ekin_GeV,
	'voltage_MV': RF_voltage_MV,
	'phase': RF_phase
}
