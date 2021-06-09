#-----------------------------------------------------------------------
# Class to output single particle co-ordinates, based on 
# Hannes Bartosik's (CERN BE-ABP-HSI) output dictionary.
# 24.07.2018: Created by Haroon Rafique, CERN BE-ABP-HSI 
# 25.05.2021: Updated by Haroon Rafique, STFC ISIS Synchrotron Group
# 
# Note use of python 2.7 nested dictionaries to store particle data
#
# This class is provided as an example of storage of PyORBIT particle 
# data for plotting of poincare sections without dumping bunch files as
# writing to file is computationally expensive
#
# Function AddNewParticle(n) adds EMPTY particle n to dictionary.
# Function Update(bunch, turn) stores data for particles at turn.
# Function PrintParticleForTurn(turn, n, filename) self-explanatory.
# Function PrintParticle(n, filename) self-explanatory.
# Function PrintAllParticles(filename) prints all particles for all
# turns for which Update() function was called.
#-----------------------------------------------------------------------

import orbit_mpi
import os

class ParticleOutputDictionary(object):
	
	def __init__(self):
		self.update_flag = 0
		self.particles = {}			# Top level dictionary : N : All data	
		self.particle_list = [] 	# Record indices of stored particles
		self.turn_list = [] 		# Record indices of stored turns
		
		self.AddNewParticle(int(0))
		
		print 'Particle_output_dictionary: Created initial particle data dictionary \'particles\'\n \tprinting for particle 0'
		print 'Format: {\'particle\': {\'turn\': {\'dE\': 0.0, \'px\': 0.0, \'py\': 0.0, \'y\': 0.0, \'x\': 0.0, \'z\': 0.0}}}'
		print self.particles
			
	def AddNewParticle(self, n):
		# Check that we haven't updated anything before adding a new particle
		if self.update_flag:
			print "Particle_output_dictionary::AddNewParticle: Particles already updated. Add particles before tracking. Aborting addition."
		else:			
			# Check that the particle isn't already present
			if n in self.particle_list:
				print "Particle_output_dictionary::AddNewParticle: Particle already added. Aborting addition."
			else:	
				# Append the new particle
				self.particles[n] = {} # First level in : N-1 : Particle Index
			
				# Add zero turn
				self.particles[n]['0'] = {}	# Second level : N-2 : Turn
			
				# Add each co-ordinate and set initial values to zero
				self.particles[n]['0']['x'] = 0.	# Third level : N-3 : x
				self.particles[n]['0']['px'] = 0.	# Third level : N-3 : px
				self.particles[n]['0']['y'] = 0.	# Third level : N-3 : y
				self.particles[n]['0']['py'] = 0.	# Third level : N-3 : py
				self.particles[n]['0']['z'] = 0.	# Third level : N-3 : z
				self.particles[n]['0']['dE'] = 0.	# Third level : N-3 : dE
				
				# Append to list of indices
				self.particle_list.append(n)
				print "Particle_output_dictionary::AddNewParticle: Particle ",n+1," added."
			
	
	def AddNParticles(self, n):
		for i in range(1, n, 1):
			self.AddNewParticle(int(i))
		
	def Update(self, bunch, turn, verbose=False):
		self.update_flag = 1
		
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
						
			for n in self.particle_list:
				# Create the turn dictionary
				self.particles[n][turn] = {}	# Second level : N-2 : Turn
			
				# self.particles[index][turn]['x'] = bunch.x(index)
				self.particles[n][turn]['x'] = bunch.x(n)
				self.particles[n][turn]['xp'] = bunch.xp(n)
				self.particles[n][turn]['y'] = bunch.y(n)
				self.particles[n][turn]['yp'] = bunch.yp(n)
				self.particles[n][turn]['z'] = bunch.z(n)
				self.particles[n][turn]['dE'] = bunch.dE(n)
				
		self.turn_list.append(turn)		
		if verbose:
			print "Particle_output_dictionary::update: Added turn %i" % (turn)
			print "Dictionary now:"
			print self.particles
				
	# Function to print 6D co-ordinates for a particle for 1 given turn
	def PrintParticleForTurn(self, turn, n, filename=None):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
			if filename is None:				
				filename = 'Particle_' + str(n) + '_turn_' + str(turn) + '.dat'
				
			# Check that the particle exists
			if n not in self.particle_list:
				print "Particle_output_dictionary::print_particle_for_turn: Particle not stored, use AddNewParticle(n) before tracking."
			else:
				# if file exists then append
				if os.path.exists(filename):
					f = open(filename,"a+")
					
										
				# if file doesn't exist create and add header
				else:
					f = open(filename,"w+")
					f.write("#ParticleID\tTurn\tx[m]\txp\ty[m]\typ\tz[m]\tdE[GeV]")
				
					f.write("\n%i\t%i\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f" % ( 	\
					n, turn, 										\
					self.particles[n][turn]['x'],		\
					self.particles[n][turn]['xp'],	\
					self.particles[n][turn]['y'],		\
					self.particles[n][turn]['yp'],	\
					self.particles[n][turn]['z'],		\
					self.particles[n][turn]['dE']		))
				f.close()

	# Function to print 6D co-ordinates for a particle for all turns
	def PrintParticle(self, n, filename=None):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
			if filename is None:				
				filename = 'Particle_' + str(n) + '.dat'
				
			# Check that the particle exists
			if n not in self.particle_list:
				print "Particle_output_dictionary::print_particle: Particle not stored, use AddNewParticle(n) before tracking."
			else:
				# if file exists then append
				if os.path.exists(filename):
					f = open(filename,"a+")
										
				# if file doesn't exist create and add header
				else:
					f = open(filename,"w+")
					f.write("#ParticleID\tTurn\tx[m]\txp\ty[m]\typ\tz[m]\tdE[GeV]")
				
				for t in self.turn_list:				
					f.write("\n%i\t%i\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f" % ( 	\
						n, t, 										\
						self.particles[n][t]['x'],		\
						self.particles[n][t]['xp'],		\
						self.particles[n][t]['y'],		\
						self.particles[n][t]['yp'],		\
						self.particles[n][t]['z'],		\
						self.particles[n][t]['dE'] 		))
				f.close()				
					
	# Function to print 6D co-ordinates for all particles for all turns
	def PrintAllParticles(self, filename=None):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
			if filename is None:				
				filename = 'Particles_all.dat'

			# if file exists then append
			if os.path.exists(filename):
				os.remove(filename)
				f = open(filename,"a")
									
			# if file doesn't exist create and add header
			else:
				f = open(filename,"w+")
				f.write("#ParticleID\tTurn\tx[m]\txp\ty[m]\typ\tz[m]\tdE[GeV]")
			
			for n in self.particle_list:	
				for t in self.turn_list:			
					f.write("\n%i\t%i\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f" % ( 	\
						n, t, 										\
						self.particles[n][t]['x'],		\
						self.particles[n][t]['xp'],		\
						self.particles[n][t]['y'],		\
						self.particles[n][t]['yp'],		\
						self.particles[n][t]['z'],		\
						self.particles[n][t]['dE']		))
			f.close()
