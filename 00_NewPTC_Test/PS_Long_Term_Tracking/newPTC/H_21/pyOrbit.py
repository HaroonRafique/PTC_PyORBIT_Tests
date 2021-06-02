import os
import sys
import time
import math
import timeit
import orbit_mpi
import numpy as np
import scipy.io as sio
from scipy.stats import moment

# Use switches in simulation_parameters.py in current folder
#-------------------------------------------------------------
from simulation_parameters import switches as s
from simulation_parameters import parameters as p

# utils
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.utils.consts import mass_proton, speed_of_light, pi

# bunch
from bunch import Bunch
from bunch import BunchTwissAnalysis, BunchTuneAnalysis
from orbit.bunch_utils import ParticleIdNumber

# diagnostics
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNodeAsChild
from orbit.diagnostics import addTeapotMomentsNodeSet, addTeapotStatLatsNodeSet

# PTC lattice
from libptc_orbit import *
from ext.ptc_orbit import PTC_Lattice
from ext.ptc_orbit import PTC_Node
from ext.ptc_orbit.ptc_orbit import setBunchParamsPTC, readAccelTablePTC, readScriptPTC
from ext.ptc_orbit.ptc_orbit import updateParamsPTC, synchronousSetPTC, synchronousAfterPTC
from ext.ptc_orbit.ptc_orbit import trackBunchThroughLatticePTC, trackBunchInRangePTC
from orbit.aperture import TeapotApertureNode

from lib.output_dictionary import *
from lib.save_bunch_as_matfile import *
from lib.suppress_stdout import suppress_STDOUT
from lib.pyOrbit_Bunch_Gather import *
from lib.pyOrbit_Tunespread_Calculator import *
from lib.pyOrbit_GenerateInitialDistribution import *
from lib.pyOrbit_PrintLatticeFunctionsFromPTC import *
from lib.pyOrbit_PTCLatticeFunctionsDictionary import *
from lib.pyOrbit_ParticleOutputDictionary import *
readScriptPTC_noSTDOUT = suppress_STDOUT(readScriptPTC)

# MPI stuff
#-----------------------------------------------------------------------
comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
print '\n\tStart PyORBIT simulation on MPI process: ', rank


# Function to check that a file isn't empty (common PTC file bug)
def is_non_zero_file(fpath):  
	print '\n\t\t\tis_non_zero_file:: Checking file ', fpath
	print '\n\t\t\tis_non_zero_file:: File exists = ', os.path.isfile(fpath)
	print '\n\t\t\tis_non_zero_file:: Size > 3 bytes = ', os.path.getsize(fpath)
	return os.path.isfile(fpath) and os.path.getsize(fpath) > 3

# Function to check and read PTC file
def CheckAndReadPTCFile(f):
	if is_non_zero_file(f): 
		readScriptPTC_noSTDOUT(f)
	else:
		print '\n\t\t\CheckAndReadPTCFile:: ERROR: PTC file ', f, ' is empty or does not exist, exiting'
		exit(0)

# Function to open TWISS_PTC_table.OUT and return fractional tunes
def GetTunesFromPTC():
	readScriptPTC_noSTDOUT('../PTC/twiss_script.ptc')
	with open('TWISS_PTC_table.OUT') as f:
		first_line = f.readline()
		Qx = (float(first_line.split()[2]))
		Qy = (float(first_line.split()[3]))
	os.system('rm TWISS_PTC_table.OUT')
	return Qx, Qy

# Create folder structure
#-----------------------------------------------------------------------
print '\n\t\tmkdir on MPI process: ', rank
from lib.mpi_helpers import mpi_mkdir_p
mpi_mkdir_p('input')
mpi_mkdir_p('bunch_output')
mpi_mkdir_p('output')
mpi_mkdir_p('lost')

# Lattice function dictionary to print closed orbit
#-----------------------------------------------------------------------
if s['Update_Twiss']:
	ptc_dictionary_file = 'input/ptc_dictionary.pkl'
	if not os.path.exists(ptc_dictionary_file):        
		PTC_Twiss = PTCLatticeFunctionsDictionary()
	else:
		with open(ptc_dictionary_file) as sid:
			PTC_Twiss = pickle.load(sid)

# Dictionary for simulation status
#-----------------------------------------------------------------------
import pickle # HAVE TO CLEAN THIS FILE BEFORE RUNNING A NEW SIMULATION
status_file = 'input/simulation_status.pkl'
if not os.path.exists(status_file):
	sts = {'turn': -1}
else:
	with open(status_file) as fid:
		sts = pickle.load(fid)
                
# Write tunes.str file for MAD-X input
#-----------------------------------------------------------------------
if not rank:
        script_name = 'PS_Lattice/tunes.str'
        if os.path.exists(script_name):  
                print 'tune file ' + script_name + ' already exists. Deleting'
                os.remove(script_name)

        f= open(script_name,"w")

        f.write('/**********************************************************************************\n')
        f.write('*                             Tunes for PTC-PyORBIT simulation\n')
        f.write('***********************************************************************************/\n')
        f.write('tune_x = 0.' + str(p['tunex'][-2:]) + ';\n')
        f.write('tune_y = 0.' + str(p['tuney'][-2:]) + ';\n')
        f.write('lattice_start = ' + str(p['transverse_plane_flag']) + ';')
        f.close()
orbit_mpi.MPI_Barrier(comm)

# Generate Lattice (MADX + PTC) - Use MPI to run on only one 'process'
#-----------------------------------------------------------------------
print '\nStart MADX on MPI process: ', rank
if not rank:
	# ~ os.system("../../../../madx-linux64_v5_02_00 < Flat_file.madx")
	os.system("../../../../madx-linux64_v5_06_01 < Flat_file.madx")
	# ~ os.system("/afs/cern.ch/eng/sl/MAD-X/pro/releases/5.02.00/madx-linux64 < Flat_file.madx")
	# ~ os.system("/afs/cern.ch/eng/sl/MAD-X/pro/releases/5.06.01//madx-linux64 < Flat_file.madx")
orbit_mpi.MPI_Barrier(comm)

# Generate PTC RF table
#-----------------------------------------------------------------------
print '\n\t\tCreate RF file on MPI process: ', rank
from lib.write_ptc_table import write_RFtable
from simulation_parameters import RFparameters as RF 
write_RFtable('input/RF_table.ptc', *[RF[k] for k in ['harmonic_factors','time','Ekin_GeV','voltage_MV','phase']])

# Initialize a Teapot-Style PTC lattice
#-----------------------------------------------------------------------
print '\n\t\tRead PTC flat file: on MPI process: ', rank
PTC_File = 'PTC-PyORBIT_flat_file.flt'
Lattice = PTC_Lattice("PS")
Lattice.readPTC(PTC_File)

print '\n\t\tRead PTC files on MPI process: ', rank
CheckAndReadPTCFile('PTC/fringe.ptc')
CheckAndReadPTCFile('PTC/time.ptc')
CheckAndReadPTCFile('PTC/ramp_cavities.ptc')

# Create a dictionary of parameters
#-----------------------------------------------------------------------
print '\n\t\tMake paramsDict on MPI process: ', rank
paramsDict = {}
paramsDict["length"]=Lattice.getLength()/Lattice.nHarm

# Add apertures
#-----------------------------------------------------------------------
print '\n\t\tAdd apertures on MPI process: ', rank
position = 0
for node in Lattice.getNodes():
	myaperturenode = TeapotApertureNode(1, 10, 10, position)
	node.addChildNode(myaperturenode, node.ENTRANCE)
	node.addChildNode(myaperturenode, node.BODY)
	node.addChildNode(myaperturenode, node.EXIT)
	position += node.getLength()

# Import a bunch and relevant parameters for it
#-----------------------------------------------------------------------
if sts['turn'] < 0:
	print '\n\t\tCreate bunch on MPI process: ', rank
	bunch = Bunch()
	setBunchParamsPTC(bunch)

	p['harmonic_number'] = Lattice.nHarm
	p['phi_s']           = 0
	p['gamma']           = bunch.getSyncParticle().gamma()
	p['beta']            = bunch.getSyncParticle().beta()
        print '\n\tBETA = ', bunch.getSyncParticle().beta()
        print '\n\tGAMMA = ', bunch.getSyncParticle().gamma()
	p['energy']          = 1e9 * bunch.mass() * bunch.getSyncParticle().gamma()
	# ~ p['bunch_length'] = p['sig_z']/speed_of_light/bunch.getSyncParticle().beta()*4
	p['bunch_length'] = p['bunch_length']
	# ~ kin_Energy = bunch.getSyncParticle().kinEnergy()

	print '\n\t\tOutput simulation_parameters on MPI process: ', rank
	for i in p:
		print '\t', i, '\t = \t', p[i]

	if s['Horizontal']:
		Particle_distribution_file = generate_initial_poincare_distributionH(p['InitialDistnSigma'], p, Lattice, output_file='input/ParticleDistribution.in', summary_file='input/ParticleDistribution_summary.txt')
	else:
		Particle_distribution_file = generate_initial_poincare_distributionV(p['InitialDistnSigma'], p, Lattice, output_file='input/ParticleDistribution.in', summary_file='input/ParticleDistribution_summary.txt')
	
	print '\nbunch_orbit_to_pyorbit on MPI process: ', rank
	bunch_orbit_to_pyorbit(paramsDict["length"], bunch.getSyncParticle().kinEnergy(), Particle_distribution_file, bunch, p['n_macroparticles'] + 1) #read in only first N_mp particles.

# Add Macrosize to bunch
#-----------------------------------------------------------------------
	bunch.addPartAttr("macrosize")
	map(lambda i: bunch.partAttrValue("macrosize", i, 0, p['macrosize']), range(bunch.getSize()))
	ParticleIdNumber().addParticleIdNumbers(bunch) # Give them unique number IDs

# Dump and save as Matfile
#-----------------------------------------------------------------------
	# ~ bunch.dumpBunch("input/mainbunch_start.dat")
	print '\n\t\tSave bunch in bunch_output/mainbunch_-000001.mat on MPI process: ', rank
	saveBunchAsMatfile(bunch, "bunch_output/mainbunch_-000001")
	print '\n\t\tSave bunch in input/mainbunch.mat on MPI process: ', rank
	saveBunchAsMatfile(bunch, "input/mainbunch")
	sts['mainbunch_file'] = "input/mainbunch"

# Create empty lost bunch
#-----------------------------------------------------------------------
	lostbunch = Bunch()
	bunch.copyEmptyBunchTo(lostbunch)
	lostbunch.addPartAttr('ParticlePhaseAttributes')
	lostbunch.addPartAttr("LostParticleAttributes")	
	saveBunchAsMatfile(lostbunch, "input/lostbunch")
	sts['lostbunch_file'] = "input/lostbunch"

# Add items to pickle parameters
#-----------------------------------------------------------------------
	sts['turns_max'] = p['turns_max']
	sts['turns_update'] = p['turns_update']
	sts['turns_print'] = p['turns_print']
	sts['circumference'] = p['circumference']

bunch = bunch_from_matfile(sts['mainbunch_file'])
lostbunch = bunch_from_matfile(sts['lostbunch_file'])
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch

# ParticleOutputDictionary for poincare sections instead of dumping bunch files
#-----------------------------------------------------------------------
particleDictionary = ParticleOutputDictionary()
particleDictionary.AddNParticles(p['n_macroparticles'])


# Add tune analysis child node
#-----------------------------------------------------
parentnode_number = 97
parentnode = Lattice.getNodes()[parentnode_number]
Twiss_at_parentnode_entrance = Lattice.getNodes()[parentnode_number-1].getParamsDict()
tunes = TeapotTuneAnalysisNode("tune_analysis")

tunes.assignTwiss(*[Twiss_at_parentnode_entrance[k] for k in ['betax','alphax','etax','etapx','betay','alphay','etay','etapy']])
tunes.assignClosedOrbit(*[Twiss_at_parentnode_entrance[k] for k in ['orbitx','orbitpx','orbity','orbitpy']])
addTeapotDiagnosticsNodeAsChild(Lattice, parentnode, tunes)

# Define twiss analysis and output dictionary
#-----------------------------------------------------------------------
print '\n\t\tbunchtwissanalysis on MPI process: ', rank
bunchtwissanalysis = BunchTwissAnalysis() #Prepare the analysis class that will look at emittances, etc.
get_dpp = lambda b, bta: np.sqrt(bta.getCorrelation(5,5)) / (b.getSyncParticle().gamma()*b.mass()*b.getSyncParticle().beta()**2)
get_bunch_length = lambda b, bta: 4 * np.sqrt(bta.getCorrelation(4,4)) / (speed_of_light*b.getSyncParticle().beta())
get_eps_z = lambda b, bta: 1e9 * 4 * pi * bta.getEmittance(2) / (speed_of_light*b.getSyncParticle().beta())

output_file = 'output/output.mat'
output = Output_dictionary()
output.addParameter('turn', lambda: turn)
#output.addParameter('epsn_x', lambda: bunchtwissanalysis.getEmittanceNormalized(0))
#output.addParameter('epsn_y', lambda: bunchtwissanalysis.getEmittanceNormalized(1))
#output.addParameter('eps_z', lambda: get_eps_z(bunch, bunchtwissanalysis))
output.addParameter('intensity', lambda: bunchtwissanalysis.getGlobalMacrosize())
output.addParameter('n_mp', lambda: bunchtwissanalysis.getGlobalCount())
#output.addParameter('D_x', lambda: bunchtwissanalysis.getDispersion(0))
#output.addParameter('D_y', lambda: bunchtwissanalysis.getDispersion(1))
#output.addParameter('bunchlength', lambda: get_bunch_length(bunch, bunchtwissanalysis))
#output.addParameter('dpp_rms', lambda: get_dpp(bunch, bunchtwissanalysis))
#output.addParameter('beta_x', lambda: bunchtwissanalysis.getBeta(0))
#output.addParameter('beta_y', lambda: bunchtwissanalysis.getBeta(1))
#output.addParameter('alpha_x', lambda: bunchtwissanalysis.getAlpha(0))
#output.addParameter('alpha_y', lambda: bunchtwissanalysis.getAlpha(1))
#output.addParameter('mean_x', lambda: bunchtwissanalysis.getAverage(0))
#output.addParameter('mean_xp', lambda: bunchtwissanalysis.getAverage(1))
#output.addParameter('mean_y', lambda: bunchtwissanalysis.getAverage(2))
#output.addParameter('mean_yp', lambda: bunchtwissanalysis.getAverage(3))
#output.addParameter('mean_z', lambda: bunchtwissanalysis.getAverage(4))
#output.addParameter('mean_dE', lambda: bunchtwissanalysis.getAverage(5))
#output.addParameter('eff_beta_x', lambda: bunchtwissanalysis.getEffectiveBeta(0))
#output.addParameter('eff_beta_y', lambda: bunchtwissanalysis.getEffectiveBeta(1))
#output.addParameter('eff_epsn_x', lambda: bunchtwissanalysis.getEffectiveEmittance(0))
#output.addParameter('eff_epsn_y', lambda: bunchtwissanalysis.getEffectiveEmittance(1))
#output.addParameter('eff_alpha_x', lambda: bunchtwissanalysis.getEffectiveAlpha(0))
#output.addParameter('eff_alpha_y', lambda: bunchtwissanalysis.getEffectiveAlpha(1))
output.addParameter('gamma', lambda: bunch.getSyncParticle().gamma())


# Pre Track Bunch Twiss Analysis & Add BunchGather outputs
#-----------------------------------------------------------------------
print '\n\t\tStart tracking on MPI process: ', rank
turn = -1
bunchtwissanalysis.analyzeBunch(bunch)

calc_moments = False # these won't work unless we use a large distribution
if calc_moments:
	moments = BunchGather(bunch, turn, p) # Calculate bunch moments and kurtosis

	# Add moments and kurtosis
	output.addParameter('sig_x', lambda: moments['Sig_x'])
	output.addParameter('sig_xp', lambda: moments['Sig_xp'])
	output.addParameter('sig_y', lambda: moments['Sig_y'])
	output.addParameter('sig_yp', lambda: moments['Sig_yp'])
	output.addParameter('sig_z', lambda: moments['Sig_z'])
	output.addParameter('sig_dE', lambda: moments['Sig_dE'])

	output.addParameter('mu_x', lambda: moments['Mu_x'])
	output.addParameter('mu_xp', lambda: moments['Mu_xp'])
	output.addParameter('mu_y', lambda: moments['Mu_y'])
	output.addParameter('mu_yp', lambda: moments['Mu_yp'])
	output.addParameter('mu_z', lambda: moments['Mu_z'])
	output.addParameter('mu_dE', lambda: moments['Mu_dE'])

	output.addParameter('min_x', lambda: moments['Min_x'])
	output.addParameter('min_xp', lambda: moments['Min_xp'])
	output.addParameter('min_y', lambda: moments['Min_y'])
	output.addParameter('min_yp', lambda: moments['Min_yp'])
	output.addParameter('min_z', lambda: moments['Min_z'])
	output.addParameter('min_dE', lambda: moments['Min_dE'])

	output.addParameter('max_x', lambda: moments['Max_x'])
	output.addParameter('max_xp', lambda: moments['Max_xp'])
	output.addParameter('max_y', lambda: moments['Max_y'])
	output.addParameter('max_yp', lambda: moments['Max_yp'])
	output.addParameter('max_z', lambda: moments['Max_z'])
	output.addParameter('max_dE', lambda: moments['Max_dE'])

	output.addParameter('kurtosis_x', lambda: moments['Kurtosis_x'])
	output.addParameter('kurtosis_xp', lambda: moments['Kurtosis_xp'])
	output.addParameter('kurtosis_y', lambda: moments['Kurtosis_y'])
	output.addParameter('kurtosis_yp', lambda: moments['Kurtosis_yp'])
	output.addParameter('kurtosis_z', lambda: moments['Kurtosis_z'])
	output.addParameter('kurtosis_dE', lambda: moments['Kurtosis_dE'])

	output.addParameter('kurtosis_x_6sig', lambda: moments['Kurtosis_x_6sig'])
	output.addParameter('kurtosis_xp_6sig', lambda: moments['Kurtosis_xp_6sig'])
	output.addParameter('kurtosis_y_6sig', lambda: moments['Kurtosis_y_6sig'])
	output.addParameter('kurtosis_yp_6sig', lambda: moments['Kurtosis_yp_6sig'])
	output.addParameter('kurtosis_z_6sig', lambda: moments['Kurtosis_z_6sig'])
	output.addParameter('kurtosis_dE_6sig', lambda: moments['Kurtosis_dE_6sig'])

start_time = time.time()
last_time = time.time()
output.addParameter('turn_time', lambda: time.strftime("%H:%M:%S"))
output.addParameter('turn_duration', lambda: (time.time() - last_time))
output.addParameter('cumulative_time', lambda: (time.time() - start_time))

# PTC_Twiss must be updated before updating output
if s['Update_Twiss']:
	PTC_Twiss.UpdatePTCTwiss(Lattice, turn)
	output.addParameter('orbit_x_min', lambda: PTC_Twiss.GetMinParameter('orbit_x', turn))
	output.addParameter('orbit_x_max', lambda: PTC_Twiss.GetMaxParameter('orbit_x', turn))
	output.addParameter('orbit_y_min', lambda: PTC_Twiss.GetMinParameter('orbit_y', turn))
	output.addParameter('orbit_y_max', lambda: PTC_Twiss.GetMaxParameter('orbit_y', turn))

output.update()

print "p['n_macroparticles'] = ", p['n_macroparticles']

for i in range(0, p['n_macroparticles'],1):
	print bunch.x(i)
particleDictionary.Update(bunch, turn, verbose=True)

if os.path.exists(output_file):
	output.import_from_matfile(output_file)

# Track
#-----------------------------------------------------------------------
print '\n\t\tStart tracking on MPI process: ', rank
start_time = time.time()
last_time = time.time()

print '\n\t\tstart time = ', start_time

for turn in range(sts['turn']+1, sts['turns_max']):
	if not rank:	last_time = time.time()

	Lattice.trackBunch(bunch, paramsDict)
	bunchtwissanalysis.analyzeBunch(bunch)  # analyze twiss and emittance
	
	if calc_moments:
		moments = BunchGather(bunch, turn, p)	# Calculate bunch moments and kurtosis
		
	if s['Update_Twiss']: 
		readScriptPTC_noSTDOUT("PTC/update-twiss.ptc") # this is needed to correclty update the twiss functions in all lattice nodes in updateParamsPTC
		updateParamsPTC(Lattice,bunch) 			# to update bunch energy and twiss functions

	if turn in sts['turns_update']:	sts['turn'] = turn

	output.update()
	particleDictionary.Update(bunch, turn)

	if turn in sts['turns_print']:
		# ~ saveBunchAsMatfile(bunch, "input/mainbunch")
		# ~ saveBunchAsMatfile(bunch, "bunch_output/mainbunch_%s"%(str(turn).zfill(6)))
		# ~ saveBunchAsMatfile(lostbunch, "lost/lostbunch_%s"%(str(turn).zfill(6)))
		output.save_to_matfile(output_file)
		if not rank:
			with open(status_file, 'w') as fid:
				pickle.dump(sts, fid)
				
particleDictionary.PrintAllParticles('Poincare.dat')
