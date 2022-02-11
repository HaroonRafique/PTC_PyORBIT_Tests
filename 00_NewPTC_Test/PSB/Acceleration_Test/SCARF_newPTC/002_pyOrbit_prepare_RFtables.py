import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from lib.write_ptc_table import write_RFtable

data_dir = 'cycle_definition/'

phases = np.load(data_dir+'LHC25PhaseProgram.npy')
ramp = np.load(data_dir+'PSBL2Ramp.npy')
voltages_raw = np.loadtxt(data_dir+'LHC1ARing3Voltages.csv', delimiter=',', skiprows=3)
FGVRFC16 = voltages_raw[:6, :2].T
FGVRFC04 = voltages_raw[:10, 4:6].T
FGVRFC02 = voltages_raw[:, 8:10].T


t  = ramp[0]/1e3
pc = ramp[1]
time_ = np.arange(t[0],t[-1],1e-3)
time_ += 0.2e-3
E0 = const.physical_constants['proton mass energy equivalent'][0]/const.e
Ekin = np.sqrt(pc**2 + E0**2) - E0
Ekin_ = np.interp(time_,t,Ekin) / 1e9


harmonics = [1,2]
t0 = time_[0]

RF_voltage_ = np.zeros((2,len(time_)))
t = FGVRFC02[0]/1e3
RF_voltage_[0,:] = np.interp(time_,t,FGVRFC02[1]*1e3) / 1e6
t = FGVRFC04[0]/1e3
RF_voltage_[1,:] = np.interp(time_,t,FGVRFC04[1]*1e3) / 1e6

RF_phase_ = np.zeros((2,len(time_)))
t = phases[0]
RF_phase_[1,:] = np.interp(time_,t,phases[1])

phase_offset = np.pi
RF_voltage_[1,:] *= 0 #-1
RF_voltage_[0,:] = 0.008 #-1
RF_phase_[0,:] += phase_offset
RF_phase_[1,:] += 2*phase_offset
# RF_phase_[1,0] -= 1.9

write_RFtable('Tables/RF_DoubleHarm_50MeV.dat', harmonics, time_-time_[0], Ekin_, RF_voltage_.T, RF_phase_.T)


if 0:
	f, axs = plt.subplots(2,1, sharex=True)
	ax = axs[0]
	ax.plot(FGVRFC02[0], FGVRFC02[1], label='C02 (voltage)')
	ax.plot(FGVRFC04[0], FGVRFC04[1], label='C04 (voltage)')
	ax.set_ylabel('voltage (kV)')
	ax.plot(phases[0]*1e3, phases[1], label='C04 (phase)')
	ax = axs[1]
	ax.set_ylabel('kinetic Energy (MeV)')
	ax.plot(time_*1e3, Ekin_*1e3)
	ax.set_xlabel('cycle time (ms)')


	f,ax = plt.subplots(1)
	ax.plot(time_,RF_voltage_[0])
	ax.plot(time_,RF_voltage_[1])

	plt.show()
