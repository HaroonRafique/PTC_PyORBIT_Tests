import math
import sys
import time
import orbit_mpi
import timeit
import numpy as np
import scipy.io as sio
import os

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

# longitudinal space charge
from orbit.space_charge.sc1d import addLongitudinalSpaceChargeNode, addLongitudinalSpaceChargeNodeAsChild, SC1D_AccNode
from spacecharge import LSpaceChargeCalc

# transverse space charge
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D
# from spacecharge import SpaceChargeCalcSliceBySlice2D

from lib.output_dictionary import *
from lib.pyOrbit_GenerateMatchedDistribution_PSB import *
from lib.save_bunch_as_matfile import *
from lib.bunch_profiles import *
from lib.suppress_stdout import suppress_STDOUT
readScriptPTC_noSTDOUT = suppress_STDOUT(readScriptPTC)

print "Start ..."
comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)


#----------------------------------------------
# Create folder sctructure
#----------------------------------------------
from lib.mpi_helpers import mpi_mkdir_p
mpi_mkdir_p('input')
mpi_mkdir_p('output')

'''
#----------------------------------------------
# Dictionary for simulation status
#----------------------------------------------
import pickle
status_file = 'input/simulation_status.pkl'
if not os.path.exists(status_file):
	sts = {'turn': -1, 'mainbunch_file': "input/mainbunch", 'lostbunch_file': "input/lostbunch"}
else: 
	with open(status_file) as fid:
		sts = pickle.load(fid)
'''



#----------------------------------------------
# Initialize a Teapot-Style PTC lattice
#----------------------------------------------
PTC_File = "PSB/madx/PTC-PyORBIT_flat_file.flt"
Lattice = PTC_Lattice("PSB")
Lattice.readPTC(PTC_File)
readScriptPTC('ptc/fringe.txt')
readScriptPTC('ptc/time.txt')
readScriptPTC('ptc/ramp_cavities.ptc')
readScriptPTC('ptc/energize_lattice.ptc')
readScriptPTC('ptc/chrom.txt')
paramsDict = {}
paramsDict["length"]=Lattice.getLength()/Lattice.nHarm


#----------------------------------------------
# Add the main bunch and lost particles bunch
#----------------------------------------------

bunch = Bunch()
setBunchParamsPTC(bunch)

from simulation_parameters import parameters as p
''
p['harmonic_number'] = Lattice.nHarm
# p['rf_voltage']      = 1e6*RF['voltage_MV'][0]
# p['phi_s']           = pi + RF['phase'][0]
p['gamma']           = bunch.getSyncParticle().gamma()
p['beta']            = bunch.getSyncParticle().beta()
p['energy']          = 1e9 * bunch.mass() * bunch.getSyncParticle().gamma()
Lattice.orbitx0 = p['orbitx0_injection']
Lattice.orbity0 = p['orbity0_injection']

directory = p['particledistribution_directory']
mpi_mkdir_p(directory)

n_macroparticles_left = np.copy(p['n_macroparticles'])
n_macroparticles_sum = 0
for turns_left in range(1,1+p['nturns_accumulation'])[::-1]:
	p['n_macroparticles'] = int(n_macroparticles_left/turns_left)
	n_macroparticles_left -= p['n_macroparticles']
	Particle_distribution_file = generate_initial_distribution_noLongPainting(p, Lattice, 
								summary_mat_file='simulation_parameters.mat',
								output_file=directory+'/OrbitL%d.dat'%(turns_left),
								summary_file=directory+'/ParticleDistribution_summary.txt', outputFormat='Orbit')
	n_macroparticles_sum += p['n_macroparticles']

print "\ntotal number of macroparticles generated: %i"%n_macroparticles_sum

# kin_Energy = 0
# bunch_orbit_to_pyorbit(paramsDict["length"], kin_Energy, Particle_distribution_file, bunch)

'''
import matplotlib.pylab as plt
dE = map(lambda i: bunch.dE(i), range(bunch.getSize()))
f, ax = plt.subplots(1)
ax.hist(dE,50)
ax.set_xlabel('dE')

z = map(lambda i: bunch.z(i), range(bunch.getSize()))
f, ax = plt.subplots(1)
ax.hist(z,50)
ax.set_xlabel('z')

plt.show()
'''



