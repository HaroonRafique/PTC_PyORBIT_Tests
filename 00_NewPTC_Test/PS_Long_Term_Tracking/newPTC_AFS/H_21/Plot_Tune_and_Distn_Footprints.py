#!/usr/bin/env python2
import os
import sys
import glob
import imageio
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')   # suppress opening of plots
import matplotlib.cm as cm
from math import log10, floor
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
from scipy.stats import moment, kurtosis

def make_directory(path):
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)  
        
class resonance_lines(object):
	
	def __init__(self, Qx_range, Qy_range, orders, periodicity):
		
		if np.std(Qx_range):
			self.Qx_min = np.min(Qx_range)
			self.Qx_max = np.max(Qx_range)
		else:
			self.Qx_min = np.floor(Qx_range)-0.05
			self.Qx_max = np.floor(Qx_range)+1.05
		if np.std(Qy_range):
			self.Qy_min = np.min(Qy_range)
			self.Qy_max = np.max(Qy_range)
		else:
			self.Qy_min = np.floor(Qy_range)-0.05
			self.Qy_max = np.floor(Qy_range)+1.05

		self.periodicity = periodicity
									
		nx, ny = [], []

		for order in np.nditer(np.array(orders)):
			t = np.array(range(-order, order+1))
			nx.extend(order - np.abs(t))
			ny.extend(t)
		nx = np.array(nx)
		ny = np.array(ny)
	
		cextr = np.array([nx*np.floor(self.Qx_min)+ny*np.floor(self.Qy_min), \
						  nx*np.ceil(self.Qx_max)+ny*np.floor(self.Qy_min), \
						  nx*np.floor(self.Qx_min)+ny*np.ceil(self.Qy_max), \
						  nx*np.ceil(self.Qx_max)+ny*np.ceil(self.Qy_max)], dtype='int')
		cmin = np.min(cextr, axis=0)
		cmax = np.max(cextr, axis=0)
		res_sum = [range(cmin[i], cmax[i]+1) for i in xrange(cextr.shape[1])]								
		self.resonance_list = zip(nx, ny, res_sum)
		
	def plot_resonance(self, figure_object = None):	
		plt.ion() # turn on interactive plotting
		if figure_object:
			fig = figure_object
			plt.figure(fig.number)
		else:
			fig = plt.figure()
		Qx_min = self.Qx_min
		Qx_max = self.Qx_max
		Qy_min = self.Qy_min
		Qy_max = self.Qy_max 
		plt.xlim(Qx_min, Qx_max)
		plt.ylim(Qy_min, Qy_max)
		plt.xlabel('Qx')
		plt.ylabel('Qy')
		for resonance in self.resonance_list:
			nx = resonance[0]
			ny = resonance[1]
			for res_sum in resonance[2]:
				if ny:
					line, = plt.plot([Qx_min, Qx_max], \
					    [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny])
				else:
					line, = plt.plot([np.float(res_sum)/nx, np.float(res_sum)/nx],[Qy_min, Qy_max])
				if ny%2:
					plt.setp(line, linestyle='--') # for skew resonances
				if res_sum%self.periodicity:
					plt.setp(line, color='b')	# non-systematic resonances
				else:
					plt.setp(line, color='r', linewidth=2.0) # systematic resonances
		plt.draw()
		return fig
        
	def plot_resonance_ax(self, axes_object, figure_object):
		plt.ioff() # turn off interactive plotting
		fig = figure_object
		ax = axes_object
		Qx_min = self.Qx_min
		Qx_max = self.Qx_max
		Qy_min = self.Qy_min
		Qy_max = self.Qy_max 
		ax.set_xlim(Qx_min, Qx_max)
		ax.set_ylim(Qy_min, Qy_max)
		ax.set_xlabel('Qx')
		ax.set_ylabel('Qy')
		for resonance in self.resonance_list:
			nx = resonance[0]
			ny = resonance[1]
			for res_sum in resonance[2]:
				if ny:
					line, = ax.plot([Qx_min, Qx_max], \
					    [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny])
				else:
					line, = ax.plot([np.float(res_sum)/nx, np.float(res_sum)/nx],[Qy_min, Qy_max])
				if ny%2:
					plt.setp(line, linestyle='--') # for skew resonances
				if res_sum%self.periodicity:
					plt.setp(line, color='b')	# non-systematic resonances
				else:
					plt.setp(line, color='r', linewidth=2.0) # systematic resonances
		fig.canvas.draw()
		ax.draw(fig.canvas.renderer, inframe=False)
		return ax
		
	def print_resonances(self):
		for resonance in self.resonance_list:
			for res_sum in resonance[2]:
				'''
				print str(resonance[0]).rjust(3), 'Qx ', ("+", "-")[resonance[1]<0], \
					  str(abs(resonance[1])).rjust(2), 'Qy = ', str(res_sum).rjust(3), \
					  '\t', ("(non-systematic)", "(systematic)")[res_sum%self.periodicity==0]
				'''
				print '%s %s%s = %s\t%s'%(str(resonance[0]).rjust(2), ("+", "-")[resonance[1]<0], \
						str(abs(resonance[1])).rjust(2), str(res_sum).rjust(4), \
						("(non-systematic)", "(systematic)")[res_sum%self.periodicity==0])
                        
def add_input_file(dd, filename, label):
	f = filename
	p = dict()
	sio.loadmat(f, mdict=p)
	dd[label] = p
	print '\tAdded output data from ', filename, '\t dictionary key: ', label
	return dd
    
def gaussian(x, A, mu, sig):
    """gaussian_3_parameter(x, A, mu, sig)"""
    return A*np.exp(-(x-mu)**2/(2*sig**2))

def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

# Plot parameters
#-----------------------------------------------------------------------
plt.close('all')

plt.rcParams['figure.figsize'] = [5.0, 4.5]
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14

plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 10

plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 5


# Location of bunch files
#-----------------------------------------------------------------------
source_dir =  'bunch_output/'
folder = source_dir

files = glob.glob(source_dir + '*.mat')
sorted(files)
files = sorted(files)[1:] # remove 0 turn

# Current directory flags
#-----------------------------------------------------------------------    
space_charge_flag = int(os.getcwd().split('/')[-1][0])
transverse_plane = os.getcwd().split('/')[-1][2]
scan_tune = os.getcwd().split('/')[-1][-2:]

main_label = ''
sc = 'SbS'
legend_label = 'Turn'
#master_bins = 128 # For mini simulation N_mp = 5E4, 64 x 64 x 32
master_bins = 512 # For full simulation N_mp = 5E5, 128 x 128 x 64

# Tune footprint
#-----------------------------------------------------------------------  
max_1d_hist = 20

min_tune = 5.80
max_tune = 6.25
q_fine = np.arange(5.5, 6.51, 0.01)
r = resonance_lines((min_tune, max_tune),(min_tune, max_tune),(1,2,3,4),10)

first_turn = False

for file in sorted(files, reverse=False):
    
    tune_tit = '(6.21, 6.245)'    

    # bunch_output/mainbunch_000001.mat
    turn = int(file.split('/')[-1].split('_')[1].split('.')[0])
    case = os.getcwd().split('/')[-1]

    if transverse_plane is 'V':
        tune_tit = '(6.21, 6.' + scan_tune + ')'
    else:
        tune_tit = '(6.' + scan_tune + ', 6.24)'
        
    print '\n\t Plotting ', case, ' scan tune =', tune_tit, ' turn = ', turn, ' tune footprint'
    
    # Load data 
    #------------------------------------------------------------------------------
    particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']

    qx = particles.ParticlePhaseAttributes[2,:]
    qy = particles.ParticlePhaseAttributes[3,:]
    qx[np.where(qx>0.5)] -= 1
    qy[np.where((qy>0.6) & (qx<0.25))] -= 1 
     
    my_cmap = plt.cm.jet
    my_cmap.set_under('w',1)

    title = str( 'Working Point = ' + tune_tit + ', Turn = ' + str(turn) )    
    
    plt.rcParams['figure.figsize'] = [6.0, 6.0]
    fig1, ax1 = plt.subplots(constrained_layout=True)
    ax1.set_title(title) 
    #r = resonance_lines((min_tune, max_tune),(min_tune, max_tune),(1,2,3,4),10)
 
    # Calculate RMS Qx and Qy
    #------------------------------------------------------------------------------
    Q_x_rms = np.sqrt(moment(6+qx,2))
    Q_y_rms = np.sqrt(moment(6+qy,2))
    Q_x_min = np.min(6+qx)
    Q_x_max = np.max(6+qx)
    Q_y_min = np.min(6+qy)
    Q_y_max = np.max(6+qy)
    Delta_q_x = Q_x_max - Q_x_min
    Delta_q_y = Q_y_max - Q_y_min    
    Delta_q_x_4sig = 4 * Q_x_rms
    Delta_q_y_4sig = 4 * Q_y_rms
    Delta_q_x_6sig = 6 * Q_x_rms
    Delta_q_y_6sig = 6 * Q_y_rms

    
    # MAIN PLOT: TUNE FOOTPRINT
    #------------------------------------------------------------------------------
    #r.plot_resonance(fig1)
    r.plot_resonance_ax(ax1, fig1)
    ax1.hist2d(6+qx, 6+qy, bins=master_bins, cmap=my_cmap, vmin=1, range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]])
    ax1.set_xlabel(r'Q$_x$')
    ax1.set_ylabel(r'Q$_y$')
    ax1.set_ylim(min_tune, max_tune)
    ax1.grid(which='both', ls=':', lw=0.5)
 
    #plt.tight_layout()
    savename = str(folder + '/Tune_Footprint_' + case + '_turn_' + str(turn) + '.png' )

    fig1.savefig(savename, dpi=500)
    plt.close(fig1)    
    
print 'Tune footprint complete'

# Vertical Phase Space
#-----------------------------------------------------------------------  
if transverse_plane is 'V':

    max_1d_hist = 20
    first_turn = True

    for file in sorted(files, reverse=True):        
    
        tune_tit = '(6.21, 6.245)'    

        # bunch_output/mainbunch_000001.mat
        turn = int(file.split('/')[-1].split('_')[1].split('.')[0])
        case = os.getcwd().split('/')[-1]

        if transverse_plane is 'V':
            tune_tit = '(6.21, 6.' + scan_tune + ')'
        else:
            tune_tit = '(6.' + scan_tune + ', 6.24)'
            
        print '\n\t Plotting ', case, ' scan tune =', tune_tit, ' turn = ', turn, ' vertical phase space'
        
        # Load data 
        #------------------------------------------------------------------------------
        particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']
        y  = particles.y *1E3
        yp = particles.yp *1E3

         
        if first_turn:
            x_max_ = np.max(y)
            x_min_ = np.min(y)
            y_min_ = np.min(yp)
            y_max_ = np.max(yp)
            if np.abs(x_max_) > np.abs(x_min_):
                x_max = round_sig(x_max_)
                x_min = round_sig(-x_max_)
            else:
                x_max = round_sig(-x_min_)
                x_min = round_sig(x_min_)
            if np.abs(y_max_) > np.abs(y_min_):
                y_max = round_sig(y_max_)
                y_min = round_sig(-y_max_)
            else:
                y_max = round_sig(-y_min_)
                y_min = round_sig(y_min_)
            first_turn = False
        
        my_cmap = plt.cm.jet
        my_cmap.set_under('w',1)

        title = str( 'Working Point = ' + tune_tit + ', Turn = ' + str(turn) )    
        
        plt.rcParams['figure.figsize'] = [6.0, 6.0]
        fig1, ax1 = plt.subplots(constrained_layout=True)
        ax1.set_title(title) 
      
        # MAIN PLOT: TUNE FOOTPRINT
        #------------------------------------------------------------------------------
        #ax1 = fig.add_subplot(gs[1:3, 0:2])
        ax1.hist2d(y, yp, bins=master_bins, cmap=my_cmap, vmin=1)
        ax1.set_xlabel('y [mm]')
        ax1.set_ylabel(r'$y^{\prime}$ [mrad]')
        ax1.set_ylim(y_min, y_max)
        ax1.set_xlim(x_min, x_max)
        ax1.grid(which='both', ls=':', lw=0.5)
     
        #plt.tight_layout()
        savename = str(folder + '/Vertical_Phase_Space_' + case + '_turn_' + str(turn) + '.png' )

        fig1.savefig(savename, dpi=500)
        plt.close(fig1)
    
    print 'Vertical phase space complete'


# Vertical Phase Space + Tune footprint
#-----------------------------------------------------------------------  
    plt.rcParams['figure.figsize'] = [8.5,4.0]

    plt.rcParams['figure.dpi'] = 200
    plt.rcParams['savefig.dpi'] = 200

    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 14

    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10

    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['lines.markersize'] = 5
    

    max_1d_hist = 20

    first_turn = True

    for file in sorted(files, reverse=True):        
    
        tune_tit = '(6.21, 6.245)'    

        # bunch_output/mainbunch_000001.mat
        turn = int(file.split('/')[-1].split('_')[1].split('.')[0])
        case = os.getcwd().split('/')[-1]

        if transverse_plane is 'V':
            tune_tit = '(6.21, 6.' + scan_tune + ')'
        else:
            tune_tit = '(6.' + scan_tune + ', 6.24)'
            
        print '\n\t Plotting ', case, ' scan tune =', tune_tit, ' turn = ', turn, ' vertical phase space and tune footprint'
        
        # Load data 
        #------------------------------------------------------------------------------
        particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']

        y  = particles.y *1E3
        yp = particles.yp *1E3
        qx = particles.ParticlePhaseAttributes[2,:]
        qy = particles.ParticlePhaseAttributes[3,:]
        qx[np.where(qx>0.5)] -= 1
        qy[np.where((qy>0.6) & (qx<0.25))] -= 1 
         
        my_cmap = plt.cm.jet
        my_cmap.set_under('w',1)

        title = str( 'Working Point = ' + tune_tit + ', Turn = ' + str(turn) )    
        
        fig1 = plt.figure(constrained_layout=True)
        gs = gridspec.GridSpec(nrows=1,ncols=3,figure=fig1,width_ratios= [1, 1E-18, 1],height_ratios=[1],wspace=0.0,hspace=0.0)
            
        if first_turn:
            x_max_ = np.max(y)
            x_min_ = np.min(y)
            y_min_ = np.min(yp)
            y_max_ = np.max(yp)
            if np.abs(x_max_) > np.abs(x_min_):
                x_max = round_sig(x_max_)
                x_min = round_sig(-x_max_)
            else:
                x_max = round_sig(-x_min_)
                x_min = round_sig(x_min_)
            if np.abs(y_max_) > np.abs(y_min_):
                y_max = round_sig(y_max_)
                y_min = round_sig(-y_max_)
            else:
                y_max = round_sig(-y_min_)
                y_min = round_sig(y_min_)
            first_turn = False    
        
        ax2 = fig1.add_subplot(gs[0, 2])
        ax2.hist2d(y, yp, bins=master_bins, cmap=my_cmap, vmin=1)
        ax2.set_xlabel('y [mm]')
        ax2.set_ylabel(r'$y^{\prime}$ [mrad]')
        ax2.set_ylim(y_min, y_max)
        ax2.set_xlim(x_min, x_max)
        ax2.grid(which='both', ls=':', lw=0.5)  
            
        # MAIN PLOT: TUNE FOOTPRINT
        #------------------------------------------------------------------------------

        ax1 = fig1.add_subplot(gs[0, 0])
        ax1.set_title(title) 
        r.plot_resonance_ax(ax1, fig1)
        #r.plot_resonance(fig1)
        ax1.hist2d(6+qx, 6+qy, bins=master_bins, cmap=my_cmap, vmin=1, range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]]) 
        ax1.set_xlabel(r'Q$_x$')
        ax1.set_ylabel(r'Q$_y$')
        ax1.set_ylim(min_tune, max_tune)
        ax1.grid(which='both', ls=':', lw=0.5)
        
        savename = str(folder + '/V_Tune_and_Phase_' + case + '_turn_' + str(turn) + '.png' )

        fig1.savefig(savename)
        plt.close(fig1)

    print 'Vertical phase space and tune footprint script complete'

# Horizontal Phase Space
#-----------------------------------------------------------------------  
else:
    max_1d_hist = 20

    first_turn = True

    for file in sorted(files, reverse=True):
        
        tune_tit = '(6.21, 6.245)'    

        # bunch_output/mainbunch_000001.mat
        turn = int(file.split('/')[-1].split('_')[1].split('.')[0])
        case = os.getcwd().split('/')[-1]

        if transverse_plane is 'V':
            tune_tit = '(6.21, 6.' + scan_tune + ')'
        else:
            tune_tit = '(6.' + scan_tune + ', 6.24)'
            
        print '\n\t Plotting ', case, ' scan tune =', tune_tit, ' turn = ', turn, ' horizontal phase space'
        
        # Load data 
        #------------------------------------------------------------------------------
        particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']
        x  = particles.x *1E3
        xp = particles.xp *1E3

         
        if first_turn:
            x_max_ = np.max(x)
            x_min_ = np.min(x)
            y_min_ = np.min(xp)
            y_max_ = np.max(xp)
            if np.abs(x_max_) > np.abs(x_min_):
                x_max = round_sig(x_max_)
                x_min = round_sig(-x_max_)
            else:
                x_max = round_sig(-x_min_)
                x_min = round_sig(x_min_)
            if np.abs(y_max_) > np.abs(y_min_):
                y_max = round_sig(y_max_)
                y_min = round_sig(-y_max_)
            else:
                y_max = round_sig(-y_min_)
                y_min = round_sig(y_min_)
            first_turn = False
        
        my_cmap = plt.cm.jet
        my_cmap.set_under('w',1)

        title = str( 'Working Point = ' + tune_tit + ', Turn = ' + str(turn) )    
        
        plt.rcParams['figure.figsize'] = [6.0, 6.0]
        fig1, ax1 = plt.subplots(constrained_layout=True)
        ax1.set_title(title) 
      
        # MAIN PLOT: TUNE FOOTPRINT
        #------------------------------------------------------------------------------
        ax1.hist2d(x, xp, bins=master_bins, cmap=my_cmap, vmin=1) 
        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel(r'$x^{\prime}$ [mrad]')
        ax1.set_ylim(y_min, y_max)
        ax1.set_xlim(x_min, x_max)
        ax1.grid(which='both', ls=':', lw=0.5)
     
        #plt.tight_layout()
        savename = str(folder + '/Horizontal_Phase_Space_' + case + '_turn_' + str(turn) + '.png' )

        fig1.savefig(savename, dpi=500)
        plt.close(fig1)
        
    print 'Horizontal phase space complete'

# Horizontal Phase Space + Tune footprint
#-----------------------------------------------------------------------  
    plt.rcParams['figure.figsize'] = [8.5,4.0]

    plt.rcParams['figure.dpi'] = 200
    plt.rcParams['savefig.dpi'] = 200

    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 14

    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10

    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['lines.markersize'] = 5

    max_1d_hist = 20

    first_turn = True

    for file in sorted(files, reverse=True):
        
        tune_tit = '(6.21, 6.245)'    

        # bunch_output/mainbunch_000001.mat
        turn = int(file.split('/')[-1].split('_')[1].split('.')[0])
        case = os.getcwd().split('/')[-1]

        if transverse_plane is 'V':
            tune_tit = '(6.21, 6.' + scan_tune + ')'
        else:
            tune_tit = '(6.' + scan_tune + ', 6.24)'
            
        print '\n\t Plotting ', case, ' scan tune =', tune_tit, ' turn = ', turn, ' horizontal phase space and tune footprint'
                
        # Load data 
        #------------------------------------------------------------------------------
        particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']
        x  = particles.x *1E3
        xp = particles.xp *1E3
        qx = particles.ParticlePhaseAttributes[2,:]
        qy = particles.ParticlePhaseAttributes[3,:]
        qx[np.where(qx>0.5)] -= 1
        qy[np.where((qy>0.6) & (qx<0.25))] -= 1 
         
        my_cmap = plt.cm.jet
        my_cmap.set_under('w',1)

        title = str( 'Working Point = ' + tune_tit + ', Turn = ' + str(turn) )    
        
        fig1 = plt.figure(constrained_layout=True)
        gs = gridspec.GridSpec(nrows=1,ncols=3,figure=fig1,width_ratios= [1, 1E-18, 1],height_ratios=[1],wspace=0.0,hspace=0.0)
            
        if first_turn:
            x_max_ = np.max(x)
            x_min_ = np.min(x)
            y_min_ = np.min(xp)
            y_max_ = np.max(xp)
            if np.abs(x_max_) > np.abs(x_min_):
                x_max = round_sig(x_max_)
                x_min = round_sig(-x_max_)
            else:
                x_max = round_sig(-x_min_)
                x_min = round_sig(x_min_)
            if np.abs(y_max_) > np.abs(y_min_):
                y_max = round_sig(y_max_)
                y_min = round_sig(-y_max_)
            else:
                y_max = round_sig(-y_min_)
                y_min = round_sig(y_min_)
            first_turn = False    
                       
        ax2 = fig1.add_subplot(gs[0, 2])
        ax2.hist2d(x, xp, bins=master_bins, cmap=my_cmap, vmin=1)
        ax2.set_xlabel('x [mm]')
        ax2.set_ylabel(r'$x^{\prime}$ [mrad]')
        ax2.set_ylim(y_min, y_max)
        ax2.set_xlim(x_min, x_max)
        ax2.grid(which='both', ls=':', lw=0.5)  
            
        # MAIN PLOT: TUNE FOOTPRINT
        #------------------------------------------------------------------------------
        ax1 = fig1.add_subplot(gs[0, 0])
        ax1.set_title(title) 
        r.plot_resonance_ax(ax1, fig1)
        #r.plot_resonance(fig1)
        ax1.hist2d(6+qx, 6+qy, bins=master_bins, cmap=my_cmap, vmin=1, range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]])
        ax1.set_xlabel(r'Q$_x$')
        ax1.set_ylabel(r'Q$_y$')
        ax1.set_ylim(min_tune, max_tune)
        ax1.grid(which='both', ls=':', lw=0.5)
        
        savename = str(folder + '/H_Tune_and_Phase_' + case + '_turn_' + str(turn) + '.png' )
        fig1.savefig(savename)
        plt.close(fig1)
        
print 'Horizontal phase space and tune footprint script complete'
  
# Plot 2x2 grid: Tune, x xp, y yp, longitudinal
#-----------------------------------------------------------------------
plt.rcParams['figure.figsize'] = [8.0, 8.0]

plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14

plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

plt.rcParams['font.size'] = 12
plt.rcParams['legend.fontsize'] = 10

plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 5

max_1d_hist = 20

min_tune = 5.80
max_tune = 6.25
q_fine = np.arange(5.5, 6.51, 0.01)

first_turn = True

for file in sorted(files, reverse=True):
        
    tune_tit = '(6.21, 6.245)'    

    # bunch_output/mainbunch_000001.mat
    turn = int(file.split('/')[-1].split('_')[1].split('.')[0])
    case = os.getcwd().split('/')[-1]

    if transverse_plane is 'V':
        tune_tit = '(6.21, 6.' + scan_tune + ')'
    else:
        tune_tit = '(6.' + scan_tune + ', 6.24)'
        
    print '\n\t Plotting ', case, ' scan tune =', tune_tit, ' turn = ', turn, ' All'
    
    # Load data 
    #------------------------------------------------------------------------------
    particles = sio.loadmat(file, squeeze_me=True,  struct_as_record=False)['particles']
    x  = particles.x*1E3     # default metres
    xp = particles.xp*1E3    # default radians
    y  = particles.y *1E3    # default metres
    yp = particles.yp *1E3   # default radians
    z  = particles.z         # default metres
    dE = particles.dE *1E3   # default GeV
    qx = particles.ParticlePhaseAttributes[2,:]
    qy = particles.ParticlePhaseAttributes[3,:]
    qx[np.where(qx>0.5)] -= 1
    qy[np.where((qy>0.6) & (qx<0.25))] -= 1 
     
    my_cmap = plt.cm.jet
    my_cmap.set_under('w',1)

    title = str( 'Working Point = ' + tune_tit + ', Turn = ' + str(turn) )    
    
    fig1 = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(nrows=2,ncols=2,figure=fig1,width_ratios= [1, 1],height_ratios=[1, 1],wspace=10.0,hspace=0.0)
    fig1.suptitle(title)
 
    # Calculate limits from largest beam (at last turn)
    if first_turn:
        x_max_1 = np.max(y)
        x_min_1 = np.min(y)
        y_min_1 = np.min(yp)
        y_max_1 = np.max(yp)
        
        if np.abs(x_max_1) > np.abs(x_min_1):
            x_max1 = round_sig(x_max_1)
            x_min1 = round_sig(-x_max_1)
        else:
            x_max1 = round_sig(-x_min_1)
            x_min1 = round_sig(x_min_1)
        if np.abs(y_max_1) > np.abs(y_min_1):
            y_max1 = round_sig(y_max_1)
            y_min1 = round_sig(-y_max_1)
        else:
            y_max1 = round_sig(-y_min_1)
            y_min1 = round_sig(y_min_1)

        x_max_2 = np.max(x)
        x_min_2 = np.min(x)
        y_min_2 = np.min(xp)
        y_max_2 = np.max(xp)
        
        if np.abs(x_max_2) > np.abs(x_min_2):
            x_max2 = round_sig(x_max_2)
            x_min2 = round_sig(-x_max_2)
        else:
            x_max2 = round_sig(-x_min_2)
            x_min2 = round_sig(x_min_2)
        if np.abs(y_max_2) > np.abs(y_min_2):
            y_max2 = round_sig(y_max_2)
            y_min2 = round_sig(-y_max_2)
        else:
            y_max2 = round_sig(-y_min_2)
            y_min2 = round_sig(y_min_2)
        
        x_max_3 = np.max(z)
        x_min_3 = np.min(z)
        y_min_3 = np.min(dE)
        y_max_3 = np.max(dE)
        
        if np.abs(x_max_3) > np.abs(x_min_3):
            x_max3 = round_sig(x_max_3)
            x_min3 = round_sig(-x_max_3)
        else:
            x_max3 = round_sig(-x_min_3)
            x_min3 = round_sig(x_min_3)
        if np.abs(y_max_3) > np.abs(y_min_3):
            y_max3 = round_sig(y_max_3)
            y_min3 = round_sig(-y_max_3)
        else:
            y_max3 = round_sig(-y_min_3)
            y_min3 = round_sig(y_min_3)
            
        first_turn = False    
    
    # x xp
    ax2 = fig1.add_subplot(gs[0, 0])
    ax2.hist2d(x, xp, bins=master_bins, cmap=my_cmap, vmin=1) 
    ax2.set_xlabel('x [mm]')
    ax2.set_ylabel(r'$x^{\prime}$ [mrad]')
    ax2.set_ylim(y_min2, y_max2)
    ax2.set_xlim(x_min2, x_max2)
    ax2.grid(which='both', ls=':', lw=0.5)  
    
    # y yp
    ax1 = fig1.add_subplot(gs[1, 1])
    ax1.hist2d(y, yp, bins=master_bins, cmap=my_cmap, vmin=1)
    ax1.set_xlabel('y [mm]')
    ax1.set_ylabel(r'$y^{\prime}$ [mrad]')
    ax1.set_ylim(y_min1, y_max1)
    ax1.set_xlim(x_min1, x_max1)
    ax1.grid(which='both', ls=':', lw=0.5)  
    
    # z dE
    ax3 = fig1.add_subplot(gs[0, 1])
    ax3.hist2d(z, dE, bins=master_bins, cmap=my_cmap, vmin=1) 
    ax3.set_xlabel('z [m]')
    ax3.set_ylabel('dE [MeV]')
    ax3.set_ylim(y_min3, y_max3)
    ax3.set_xlim(x_min3, x_max3)
    ax3.grid(which='both', ls=':', lw=0.5)  
    
    #ax1 = fig.add_subplot(gs[1:3, 0:2])
        
    # MAIN PLOT: TUNE FOOTPRINT
    #------------------------------------------------------------------------------
    ax4 = fig1.add_subplot(gs[1, 0])
    #r.plot_resonance(fig1)
    r.plot_resonance_ax(ax4, fig1)
    ax4.hist2d(6+qx, 6+qy, bins=master_bins, cmap=my_cmap, vmin=1, range=[[r.Qx_min, r.Qx_max], [r.Qy_min, r.Qy_max]]) 
    ax4.set_xlabel(r'Q$_x$')
    ax4.set_ylabel(r'Q$_y$')
    ax4.set_ylim(min_tune, max_tune)
    ax4.grid(which='both', ls=':', lw=0.5)
    
    savename = str(folder + '/All_' + case + '_turn_' + str(turn) + '.png' )

    fig1.savefig(savename)
    plt.close(fig1)
    
print 'Plotting script complete'
