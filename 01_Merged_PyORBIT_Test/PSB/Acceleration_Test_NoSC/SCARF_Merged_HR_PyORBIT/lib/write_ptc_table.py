import numpy as np
import orbit_mpi
from scipy.misc import factorial

def write_RFtable(filename, harmonic_factors, time, E_kin, RF_voltage, RF_phase):
    comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
    rank = orbit_mpi.MPI_Comm_rank(comm)
    if not rank:
        n_lines = len(time)
        n_harmonics = len(harmonic_factors)
        arr = np.vstack((time,E_kin, np.dstack((RF_voltage,RF_phase)).flatten().reshape(n_lines, 2*n_harmonics).T)).T    
        with open(filename, 'w') as fid:
            fid.write('%d  1  1  0  %d\n'%(n_lines, n_harmonics))
            fid.write('  '.join(map(lambda i: '%d'%i, harmonic_factors))+'\n')
            for j in xrange(n_lines):
                fid.write('\t'.join(map(lambda i: '%1.8f'%i, arr[j, :]))+'\n')
    orbit_mpi.MPI_Barrier(comm)

def write_PTCtable(filename, multipole_orders, time, normal_components, skew_components):
    comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
    rank = orbit_mpi.MPI_Comm_rank(comm)
    if not rank:
        multipole_orders = np.atleast_1d(multipole_orders)
        factors = 1./factorial(multipole_orders-1) # the factorial factor is needed to be consistent with MADX
        normal_components = (factors.T * np.atleast_2d(normal_components).T)
        skew_components   = (factors.T * np.atleast_2d(skew_components).T)
        arr = np.empty((normal_components.shape[0], 1+normal_components.shape[1]*2), dtype=normal_components.dtype)
        arr[:,0] = time
        arr[:,1::2] = normal_components
        arr[:,2::2] = skew_components
        n_lines = len(time)
        n_multipoles = len(multipole_orders) # number of multipole orders to be changed
        with open(filename, 'w') as fid:
            fid.write('%d  1  %d\n'%(n_lines, n_multipoles))
            fid.write(' '.join(map(lambda i: '%d'%i, multipole_orders)) + '\n')
            for j in xrange(n_lines):
                fid.write('\t'.join(map(lambda i: '%+1.11f'%i, arr[j, :]))+'\n')
    orbit_mpi.MPI_Barrier(comm)
    