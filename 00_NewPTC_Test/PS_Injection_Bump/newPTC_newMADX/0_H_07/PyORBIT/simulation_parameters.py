import os
import numpy as np

# If we are using the folder naming convention:
# Space charge flag - horizontal/vertical scan - scan point
# e.g. 0_H_07

space_charge_flag = int(os.getcwd().split('/')[-2][0])
transverse_plane = os.getcwd().split('/')[-2][2]
scan_tune = os.getcwd().split('/')[-2][-2:]

# ~ space_charge_flag = 0
# ~ transverse_plane = 'H'
# ~ scan_tune = '07'

print 'simulation_parameters: space charge = ', space_charge_flag
print 'simulation_parameters: transverse_plane = ', transverse_plane
print 'simulation_parameters: scan_tune = ', scan_tune

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
        print 'For example: 0_H_07 will launch a no space charge simulation for the horizontal tune of 6.07'
        exit(0)
                
parameters['n_macroparticles']			= int(5E4)

# Make sure to fix the initial distribution at the nominal working point (6.21, 6.24)
parameters['tomo_file'] = '../../../00_Longitudinal_Distribution/PyORBIT_Tomo_file_MD4224_HB.mat'
parameters['LatticeFile'] = '../PTC_Twiss/1.ptc'
parameters['input_distn']=''

parameters['gamma']			= 2.49253731343
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

parameters['turns_max'] = int(2200)

# Define how often we dump bunch output files
#-----------------------------------------------------------------------
tu1 = range(-1, parameters['turns_max'], 100) # every 100 turns
tu2 = range(1, 100)              # every turn for the first 100 turns
tu = tu2 + tu1

parameters['turns_max'] = int(2200)
tu1                     = range(-1, parameters['turns_max'], 50)
tu2                     = range(50, 100, 10) 
tu3                     = range(1, 50) # every turn for the first 50 turns
tu                      = tu2 + tu1 + tu3 
tu.append(874) # Wire Scanner at ctime = 172 s
tu.append(2185)# Wire Scanner at ctime = 175 s
# ~ tu.append(6556)# Wire Scanner at ctime = 185 s
parameters['turns_print'] = sorted(tu)
parameters['turns_update'] = sorted(tu)

# Simulation switches
#-----------------------------------------------------------------------
switches = {
	'InjectionBump': True, # Apply the injection bump in MAD-X and load corresponding PTC tables
	'CreateDistn': False,  # Load from file to fix initial distribution to the incoming PSB beam
	'Update_Twiss':	True,  # Perform PTC twiss and dump each turn - needed to output tune changes
	'GridSizeX': 64,
	'GridSizeY': 64,
	'GridSizeZ': 32
}

if space_charge_flag:
        switches['Space_Charge'] = True        
        if transverse_plane == 'H': 
                parameters['input_distn'] = '../../../00_Longitudinal_Distribution/500000/1_H_21/bunch_output/mainbunch_-000001.mat'
                # ~ parameters['input_distn'] = '../../../07_Dispersion_Mismatch/Initial_Distributions/' + str(parameters['n_macroparticles']) + '/1_H_21/bunch_output/mainbunch_-000001.mat'
        else:
                parameters['input_distn'] = '../../../00_Longitudinal_Distribution/500000/1_V_24/bunch_output/mainbunch_-000001.mat'
                # ~ parameters['input_distn'] = '../../../07_Dispersion_Mismatch/Initial_Distributions/' + str(parameters['n_macroparticles']) + '/1_V_24/bunch_output/mainbunch_-000001.mat'
else:
        switches['Space_Charge'] = False
        if transverse_plane == 'H': 
                parameters['input_distn'] = '../../../00_Longitudinal_Distribution/500000/0_H_21/bunch_output/mainbunch_-000001.mat'
                # ~ parameters['input_distn'] = '../../../07_Dispersion_Mismatch/Initial_Distributions/' + str(parameters['n_macroparticles']) + '/0_H_21/bunch_output/mainbunch_-000001.mat'
        else:
                parameters['input_distn'] = '../../../00_Longitudinal_Distribution/500000/0_V_24/bunch_output/mainbunch_-000001.mat'
                # ~ parameters['input_distn'] = '../../../07_Dispersion_Mismatch/Initial_Distributions/' + str(parameters['n_macroparticles']) + '/0_V_24/bunch_output/mainbunch_-000001.mat'

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
