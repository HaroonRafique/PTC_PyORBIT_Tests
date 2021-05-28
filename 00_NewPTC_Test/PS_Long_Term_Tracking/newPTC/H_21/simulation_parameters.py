import os
import numpy as np

# Folder naming convention:
# Space charge flag - horizontal/vertical scan - scan point
# e.g. 0_H_07

print 'simulation_parameters: space charge not included in this example of long term tracking'
transverse_plane = os.getcwd().split('/')[-1][0]
print 'simulation_parameters: transverse_plane = ', transverse_plane
scan_tune = os.getcwd().split('/')[-1][-2:]

parameters = {}

# Nominal working point for BCMS PS beam
parameters['tunex']			= '6.21'
parameters['tuney']			= '6.24'

if transverse_plane == 'H':
        parameters['transverse_plane_flag'] = 1
        parameters['lattice_start'] 	= 'PR.BWSH65'
        parameters['tunex']             = '6.' + str(scan_tune)
elif transverse_plane == 'V':
        parameters['transverse_plane_flag'] = 0
        parameters['lattice_start'] 	= 'PR.BWSV64'
        parameters['tuney']             = '6.' + str(scan_tune)
else:
        print 'simulation_parameters: transverse plane selection not recognised, please use either H or V in the folder name'
        print 'For example: H_07 will launch a simulation for the horizontal tune of 6.07'
        exit(0)
                
parameters['turns_max'] 		= int(5E4)
parameters['n_macroparticles']	= int(21)
parameters['InitialDistnSigma'] = 10.	# Poincare distribution limits in units of sigma

parameters['gamma']				= 2.49253731343
parameters['intensity']			= 72.5E+10
parameters['bunch_length']		= 140e-9
parameters['blength']			= 140e-9
parameters['epsn_x']			= 1E-6
parameters['epsn_y']			= 1E-6
parameters['dpp_rms']			= 8.7e-04
parameters['LongitudinalJohoParameter'] = 1.2
parameters['LongitudinalCut'] 	        = 2.4
parameters['TransverseCut']		= 5
parameters['rf_voltage']		= 0.0212
parameters['circumference']		= 2*np.pi*100
parameters['phi_s']			= 0
parameters['macrosize']			= parameters['intensity']/float(parameters['n_macroparticles'])

# PS Injection 1.4 GeV
parameters['gamma'] 	= 2.49253731343
parameters['beta'] 	= np.sqrt(parameters['gamma']**2-1)/parameters['gamma']
print 'beta = ', parameters['beta'] 
c 			= 299792458
parameters['sig_z'] 	= (parameters['beta'] * c * parameters['blength'])/4.

# Define how often we dump bunch output files
#-----------------------------------------------------------------------
n_turns = 1
parameters['turns_print'] = range(-1, parameters['turns_max'], n_turns) # every n_turns
parameters['turns_update'] = range(-1, parameters['turns_max'], n_turns) # every n_turns

# Simulation switches
#-----------------------------------------------------------------------
switches = {
	'Horizontal': True,   # Initial Poincare distribution in horizontal (true) / vertical (false)
	'Update_Twiss':	False   # Perform PTC twiss and dump each turn - needed to output tune changes
}

# PTC RF Table Parameters
#-----------------------------------------------------------------------
harmonic_factors = [1] # this times the base harmonic defines the RF harmonics (for SPS = 4620, PS 10MHz 7, 8, or 9)
time = np.array([0,1,2])
ones = np.ones_like(time)
Ekin_GeV = 1.4*ones
RF_voltage_MV = np.array([parameters['rf_voltage']*ones]).T # in MV
RF_phase = np.array([np.pi*ones]).T

RFparameters = {
	'harmonic_factors': harmonic_factors,
	'time': time,
	'Ekin_GeV': Ekin_GeV,
	'voltage_MV': RF_voltage_MV,
	'phase': RF_phase
}
