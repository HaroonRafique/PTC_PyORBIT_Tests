# Plots all available plottable data from output.mat as individual files
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy.io as sio 
import matplotlib.cm as cm

plt.rcParams['figure.figsize'] = [8.0, 4.0]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 'medium'
plt.rcParams['lines.linewidth'] = 1


'''
add_input_file:
dd: dictionary of particle data dictionaries. key = user defined label
filename: input file name (with relative path from this file)
label: file label e.g. 'case 1', 'case 2', ...
'''
def add_input_file(dd, filename, label):
	f = filename
	p = dict()
	sio.loadmat(f, mdict=p)
	dd[label] = p	
	print '\tAdded output data from ', filename, '\t dictionary key: ', label
	return dd

# Function takes the data, creates the correct timing for each value in
# seconds, then appends the required number of intervals with zeroes
def Create_Timing(ramp_stop_time, simulation_stop_time):

	time = np.linspace(0, ramp_stop_time, d_len)

	# Calculate interval to complete data until the end of the simulation
	interval = ramp_stop_time / (d_len-1)

	# Append zeroes (bump is closed)
	steps = int((simulation_stop_time - ramp_stop_time)	/ interval) + 5

	for i in xrange(steps):
		if i is 0:
			pass
		else:
			time = np.append(time, (ramp_stop_time + (i*interval)))

	result = np.column_stack([time])

	return result

'''
------------------------------------------------------------------------
						Open files and read data
------------------------------------------------------------------------
'''
	
# Create dd dictionary
dd = dict()
dd = add_input_file(dd, 'output/output.mat', 'PyORBIT')
print 'Final data dictionary keys: ', dd.keys()

sc = 'Slice-by-Slice'
main_label = '08'
main_label2 = main_label + '_zoom'
scaled_label = main_label + '_scaled'
legend_label = r'$Q_y$'
turn_tot = None
zoom_turns = 15
turns = [0, 1, 10, 100, 199, 874, 2185]

madx_qx = []
madx_qy = []
ptc_qx = []
ptc_qy = []


# Open file
infile = '../MADX_Tables/BSEXT_Out.tfs'
fin = open(infile,'r').readlines()[8:]

# Save s, x
for l in fin:
	madx_qx.append(float(l.split()[9])-6.)
	madx_qy.append(float(l.split()[10])-6.)
	ptc_qx.append(float(l.split()[11]))
	ptc_qy.append(float(l.split()[12]))

size = len(ptc_qx)
turns = np.linspace(0, 500, size)

'''
------------------------------------------------------------------------
								Plot
------------------------------------------------------------------------
'''
fig1 = plt.figure(facecolor='w', edgecolor='k')
ax1 = fig1.add_subplot(111)

colors = cm.rainbow(np.linspace(0, 1, len(dd.keys())+2 ))
c_it = int(0)

ax1.set_title('Vertical Tune');

for key, value in sorted(dd.iteritems()):
	ax1.plot(dd[key]['turn'][0], dd[key]['Qy'][0], label=key, color=colors[c_it]);
	c_it = c_it + 1

ax1.plot(turns, madx_qy, label='MAD-X', color=colors[c_it]);
c_it = c_it + 1

ax1.plot(turns, ptc_qy, label='PTC', color=colors[c_it]);
c_it = c_it + 1

ax1.set_ylabel(r'$Q_y$ [-]');
ax1.set_xlabel('Turn [-]');
ax1.grid(True);
ax1.legend();

figname = 'MADX_PyORBIT_Vertical_Tunes.png'
fig1.savefig(figname);


fig1 = plt.figure(facecolor='w', edgecolor='k')
ax1 = fig1.add_subplot(111)

colors = cm.rainbow(np.linspace(0, 1, len(dd.keys())+2 ))
c_it = int(0)

ax1.set_title('Horizontal Tune');

for key, value in sorted(dd.iteritems()):
	ax1.plot(dd[key]['turn'][0], dd[key]['Qx'][0], label=key, color=colors[c_it]);
	c_it = c_it + 1

ax1.plot(turns, madx_qx, label='MAD-X', color=colors[c_it]);
c_it = c_it + 1

ax1.plot(turns, ptc_qx, label='PTC', color=colors[c_it]);
c_it = c_it + 1

ax1.set_ylabel(r'$Q_x$ [-]');
ax1.set_xlabel('Turn [-]');
ax1.grid(True);
ax1.legend();

figname = 'MADX_PyORBIT_Horizontal_Tunes.png'
fig1.savefig(figname);

