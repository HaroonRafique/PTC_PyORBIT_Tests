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
from spacecharge import SpaceChargeCalc2p5D, Boundary2D

# apertures 
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject
from orbit.aperture import addTeapotApertureNode
from orbit.aperture import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from orbit.aperture import addCircleApertureSet, addEllipseApertureSet, addRectangleApertureSet

# collimator
from orbit.collimation import TeapotCollimatorNode
from orbit.collimation import addTeapotCollimatorNode
from collimator import Collimator

# dictionary
from lib.output_dictionary import *
from lib.save_bunch_as_matfile import *
from lib.bunch_profiles import *
from lib.suppress_stdout import suppress_STDOUT
readScriptPTC_noSTDOUT = suppress_STDOUT(readScriptPTC)

print "Start ..."
comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)

#----------------------------------------------
# Create folder structure
#----------------------------------------------`
from lib.mpi_helpers import mpi_mkdir_p
mpi_mkdir_p('input')
mpi_mkdir_p('output')


#----------------------------------------------
# Dictionary for simulation status (for resume)
#----------------------------------------------
import pickle
status_file = 'input/simulation_status.pkl'
if not os.path.exists(status_file):
	sts = {'turn': -1, 'mainbunch_file': "input/mainbunch", 'lostbunch_file': "input/lostbunch"}
else: 
	with open(status_file) as fid:
		sts = pickle.load(fid)


#----------------------------------------------
# Simulation Parameters
#----------------------------------------------
from simulation_parameters import parameters as p
sts['turns_max'] = p['turns_max']
sts['turns_print'] = p['turns_print']
sts['turns_injection'] = p['turns_injection']


#----------------------------------------------
# Initialize a Teapot-Style PTC lattice
#----------------------------------------------
PTC_File='PSB/madx/PTC-PyORBIT_flat_file.flt' 
Lattice = PTC_Lattice("PSB")
Lattice.readPTC(PTC_File)
readScriptPTC('ptc/fringe.txt')
readScriptPTC('ptc/time.txt')
readScriptPTC('ptc/chrom.txt')
readScriptPTC('ptc/ramp_cavities.ptc')
if sts['turn'] >= 0:
	readScriptPTC('ptc/read_FINAL_SETTINGS.ptc')
readScriptPTC('ptc/energize_lattice.ptc')
readScriptPTC('ptc/twiss_script.ptc')


paramsDict = {}
paramsDict["length"]=Lattice.getLength()/Lattice.nHarm
print '\nLattice parameters ...'
print '  circumference: \t', Lattice.getLength(), 'm'
print '  alphax0: \t\t', Lattice.alphax0
print '  betax0: \t\t', Lattice.betax0, 'm'
print '  alphay0: \t\t', Lattice.alphay0
print '  betay0: \t\t', Lattice.betay0, 'm'
print '  Dx0: \t\t\t', Lattice.etax0, 'm'
print '  Dpx0: \t\t', Lattice.etapx0, 'm'
print '  harm. number: \t', Lattice.nHarm
print '  nodes: \t\t', Lattice.nNodes

#----------------------------------------------
# Add apertures
#----------------------------------------------
position = 0
pos_start = 0
pos_stop  = Lattice.getLength()
n=0
for node in Lattice.getNodes():
	myaperturenode = TeapotApertureNode(1, 10, 18, position)
	node.addChildNode(myaperturenode, node.ENTRANCE)
	node.addChildNode(myaperturenode, node.BODY)
	node.addChildNode(myaperturenode, node.EXIT)
	position += node.getLength()
	n += 1

#-----------------------------------------------------
# Add tune analysis child node
#-----------------------------------------------------
parentnode_number = 40
parentnode = Lattice.getNodes()[parentnode_number]
Twiss_at_parentnode_entrance = Lattice.getNodes()[parentnode_number-1].getParamsDict
tunes = TeapotTuneAnalysisNode("tune_analysis")
tunes.assignTwiss(*[Twiss_at_parentnode_entrance()[k] for k in ['betax','alphax','etax','etapx','betay','alphay','etay','etapy']])
tunes.assignClosedOrbit(*[Twiss_at_parentnode_entrance()[k] for k in ['orbitx','orbitpx','orbity','orbitpy']])
addTeapotDiagnosticsNodeAsChild(Lattice, parentnode, tunes)


#----------------------------------------------
# Add the main bunch and lost particles bunch
#----------------------------------------------
print '\nAdding main bunch ...'

macrosize = p['macrosize']
bunch = Bunch()
setBunchParamsPTC(bunch)
kin_Energy = bunch.getSyncParticle().kinEnergy()
print '  Momentum: ', bunch.getSyncParticle().momentum(), 'GeV'
print '  Ekin:     ', bunch.getSyncParticle().kinEnergy(), 'GeV'
print '  Gamma:    ',bunch.getSyncParticle().gamma() 
print '  Beta:     ', bunch.getSyncParticle().beta()
print '  Charge:   ', bunch.charge(), 'e'
print '  Mass:     ', bunch.mass(), 'GeV'

ParticleIdNumber().addParticleIdNumbers(bunch) # Give particles unique number ids
bunch.addPartAttr("macrosize")

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes")

paramsDict["bunch"]= bunch
paramsDict["lostbunch"] = lostbunch

# read the particles from main bunch if the simulation resumes
if sts['turn'] >= 0:
	bunch = bunch_from_matfile(sts['mainbunch_file'])
	lostbunch = bunch_from_matfile(sts['lostbunch_file'])
	lostbunch.addPartAttr("LostParticleAttributes")

'''
#----------------------------------------------------
# Add transverse potential space charge node with
# rectangular boundary
#----------------------------------------------------
print '\nAdding longitudinal space charge ...'
b_a = 1.5
nMacrosMin = 32
useSpaceCharge = 1
nBins= 64     # Number of longitudinal slices in the 1D space charge solver
position = 1   # The location in the lattice. Can be any empty place
length = Lattice.getLength() 
sc1Dnode = SC1D_AccNode(b_a,length, nMacrosMin, useSpaceCharge, nBins)
nodes = Lattice.getNodes()
AccNode = nodes[1]
addLongitudinalSpaceChargeNodeAsChild(Lattice, AccNode, sc1Dnode)
''
print '\nAdding transverse space charge nodes ...'
sizeX = 128
sizeY = 128
sizeZ = 64  # Number of longitudinal slies in the 2.5D solver
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)
sc_path_length_min = 0.00000001

# Add the space charge solver to the lattice as child nodes
sc_nodes = scLatticeModifications.setSC2p5DAccNodes(Lattice, sc_path_length_min, calc2p5d)
print '  Installed %i space charge nodes'%(len(sc_nodes))
'''
#----------------------------------------------------
# Define twiss analysis and output dictionary
#----------------------------------------------------
bunchtwissanalysis = BunchTwissAnalysis() #Prepare the analysis class that will look at emittances, etc.
lostbunchtwissanalysis = BunchTwissAnalysis() #Prepare the analysis class that will look at emittances, etc.
get_dpp = lambda b, bta: np.sqrt(bta.getCorrelation(5,5)) / (b.getSyncParticle().gamma()*b.mass()*b.getSyncParticle().beta()**2)
get_bunch_length = lambda b, bta: 4 * np.sqrt(bta.getCorrelation(4,4)) / (speed_of_light*b.getSyncParticle().beta())
get_eps_z = lambda b, bta: 1e9 * 4 * pi * bta.getEmittance(2) / (speed_of_light*b.getSyncParticle().beta())
def get_BunchingFactor(grid1D,bucketlength):
	maxValue = max(map(grid1D.getValueOnGrid, range(grid1D.getSizeZ())))
	sumValue = grid1D.getSum()
	return (sumValue/bucketlength)/(maxValue/grid1D.getStepZ())

output_file = 'output/output.mat'
output = Output_dictionary()
output.addParameter('turn', lambda: turn)
output.addParameter('intensity', lambda: bunchtwissanalysis.getGlobalMacrosize())
output.addParameter('n_mp', lambda: bunchtwissanalysis.getGlobalCount())
output.addParameter('intensity_lost', lambda: macrosize*lostbunchtwissanalysis.getGlobalCount())
# output.addParameter('bunching_factor', lambda: get_BunchingFactor(calc2p5d.getLongGrid(),paramsDict["length"]))
output.addParameter('gamma', lambda: bunch.getSyncParticle().gamma())
output.addParameter('beta', lambda: bunch.getSyncParticle().beta())
output.addParameter('mean_x', lambda: bunchtwissanalysis.getAverage(0))
output.addParameter('mean_xp', lambda: bunchtwissanalysis.getAverage(1))
output.addParameter('mean_y', lambda: bunchtwissanalysis.getAverage(2))
output.addParameter('mean_yp', lambda: bunchtwissanalysis.getAverage(3))
output.addParameter('mean_z', lambda: bunchtwissanalysis.getAverage(4))
output.addParameter('mean_dE', lambda: bunchtwissanalysis.getAverage(5))
output.addParameter('epsn_x', lambda: bunchtwissanalysis.getEmittanceNormalized(0))
output.addParameter('epsn_y', lambda: bunchtwissanalysis.getEmittanceNormalized(1))
output.addParameter('eps_z', lambda: get_eps_z(bunch, bunchtwissanalysis))
output.addParameter('bunchlength', lambda: get_bunch_length(bunch, bunchtwissanalysis))
output.addParameter('dpp_rms', lambda: get_dpp(bunch, bunchtwissanalysis))
output.addParameter('dE_rms', lambda: np.sqrt(bunchtwissanalysis.getCorrelation(5,5)))
output.addParameter('bunchlength_rms', lambda: np.sqrt(bunchtwissanalysis.getCorrelation(4,4)))
output.addParameter('CO_x', lambda: Lattice.orbitx0)
output.addParameter('CO_y', lambda: Lattice.orbity0)
output.addParameter('wall_time', time.time)
if os.path.exists(output_file):
	output.import_from_matfile(output_file)

#----------------------------------------------------
# Bunch profies
#----------------------------------------------------
bunchProfiles = Bunch_Profile(50,50,50)
twiss_dict = {'betax': Lattice.betax0, 'betay': Lattice.betay0, 'etax': Lattice.etax0, 'etay': Lattice.etay0}


#----------------------------------------------------
# Tracking
#----------------------------------------------------
print '\n\n now start tracking ...'

for turn in range(sts['turn']+1, sts['turns_max']):

	if turn in sts['turns_injection']:
		Particle_distribution_file = p['particledistribution_directory']+'/OrbitL%d.dat'%(turn+1)	# final distribution with the correct angle
		kin_Energy = bunch.getSyncParticle().kinEnergy()
		bunch_orbit_to_pyorbit(paramsDict["length"], kin_Energy, Particle_distribution_file, bunch) #read in N_mp particles. 
		for i in range(bunch.getSize()):
			bunch.partAttrValue("macrosize", i, 0, macrosize)  #main bunch has finite macrosize for space charge
	
	# keep particles within circumference
	z_lim = -25*np.pi #-65 #0#-25*pi*2 #-25*pi #-65
	for i in xrange(bunch.getSize()):
		bunch.z(i, ((bunch.z(i)-z_lim)%paramsDict["length"])+z_lim)

	Lattice.trackBunch(bunch, paramsDict)
	bunchtwissanalysis.analyzeBunch(bunch)  # analyze twiss and emittance	
	lostbunchtwissanalysis.analyzeBunch(lostbunch)
	readScriptPTC_noSTDOUT("ptc/update-twiss.ptc") # this is needed to correclty update the twiss functions in all lattice nodes in updateParamsPTC
	updateParamsPTC(Lattice,bunch) # to update bunch energy and twiss functions
	tunes.assignTwiss(*[Twiss_at_parentnode_entrance()[k] for k in ['betax','alphax','etax','etapx','betay','alphay','etay','etapy']])
	tunes.assignClosedOrbit(*[Twiss_at_parentnode_entrance()[k] for k in ['orbitx','orbitpx','orbity','orbitpy']])

	output.update()
	sts['turn'] = turn

	# subtract circumference each turn in order to reconstruct the turn number from loss position
	map(lambda i: lostbunch.partAttrValue("LostParticleAttributes", i, 0, 
					  lostbunch.partAttrValue("LostParticleAttributes", i, 0)-paramsDict["length"]), xrange(lostbunch.getSize()))

	if turn in sts['turns_print']:
		output.save_to_matfile(output_file)
		saveBunchAsMatfile(bunch, sts['mainbunch_file'])
		saveBunchAsMatfile(lostbunch, sts['lostbunch_file'])
		bunchProfiles.bin_bunch(bunch, twiss_dict)
		bunchProfiles.save_bunchprofile_to_matfile('output/bunch_profile_%s.mat'%(str(turn).zfill(6)))
		if not rank:
			readScriptPTC_noSTDOUT('ptc/write_FINAL_SETTINGS.ptc')
			os.system('cp %s %s'%(sts['mainbunch_file']+'.mat', "output/mainbunch_%s.mat"%(str(turn).zfill(6))))
			os.system('cp %s %s'%(sts['lostbunch_file']+'.mat', "output/lostbunch_%s.mat"%(str(turn).zfill(6))))
			with open(status_file, 'w') as fid:
				pickle.dump(sts, fid)


# make sure simulation terminates properly
orbit_mpi.MPI_Barrier(comm)
