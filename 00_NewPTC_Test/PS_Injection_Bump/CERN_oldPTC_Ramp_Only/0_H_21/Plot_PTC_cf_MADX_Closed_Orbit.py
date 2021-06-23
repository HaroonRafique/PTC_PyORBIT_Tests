# Script to read multiple PTC twiss files and plot the closed orbit
# Expects a PTC twiss created with the following command:
# select, flag=ptc_twiss, column=name, s, betx, px, bety, py, disp3, disp3p, disp1, disp1p, x, y;
# 26.08.19 Haroon Rafique CERN BE-ABP-HSI 

import matplotlib
matplotlib.use('Agg')   # suppress opening of plots
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import numpy as np
import os
import scipy.io as sio

plt.rcParams['figure.figsize'] = [5.0, 4.5]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14

plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 8

plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 5

print ('\n\PTC MAD-X Injection Bump Comparison Plotting Script Started\n')


ptc_extensions = ('.ptc')		# All outputs are .ptc files
ptc_iterators = []				# Integers (turn) used to iterate over files

# Search from current directory
print ('\nFind PTC Twiss files\n')
for subdir, dirs, files in os.walk('PTC_Twiss'):
	# Iterate over all files
    for file in files:
		# Find files with required extension
		ext = os.path.splitext(file)[-1].lower()
		if ext in ptc_extensions:
			# Print full file path
			# ~ print (os.path.join(subdir, file))		# full path to file
			fileno = int(file.split('.')[0])	# use turn number as a key
			ptc_iterators.append(fileno)

ptc_data = dict()

# iterate over files (turns)
print ('\nRead s, x ptc_data from files\n')
ptc_last_s = 0

#* NAME S BETX BETY ALFX ALFY DISP1 DISP3 X PX Y PY L DISP2 DISP4 MU1 MU2 
for turn in sorted(ptc_iterators):
	s = []
	x = []

	# Open file
	infile = 'PTC_Twiss/' + str(turn) + '.ptc'
	fin = open(infile,'r').readlines()[90:]

	# Save s, x
	for l in fin:
		if ptc_last_s == float(l.split()[1]):
			pass
		else:
			ptc_last_s =float(l.split()[1])
			s.append(float(l.split()[1]))
			x.append(float(l.split()[6])*1E3)

	# Add to dictionary as dict[turn] = (s, x)
	ptc_data[turn] = [s, x]
	ptc_last_s = 0

# Indicate location of BSW's on zoom plots
BSW40 = 0.
BSW42 = 0.
BSW43 = 0.
BSW44 = 0.

# Find BSW40, 42, 43, 44 locations
fin0=open('optimised_flat_file.tfs','r').readlines()
for l in fin0:
        if 'BSW40' in l:
                BSW40 = float(l.split()[1])
        elif 'BSW42' in l:
                BSW42 = float(l.split()[1])
        elif 'BSW43' in l:
                BSW43 = float(l.split()[1])
        elif 'BSW44' in l:
                BSW44 = float(l.split()[1])
                
print ('\n\tBSW40 found at s = ', BSW40)
print ('\n\tBSW42 found at s = ', BSW42)
print ('\n\tBSW43 found at s = ', BSW43)
print ('\n\tBSW44 found at s = ', BSW44)
                
madx_extensions = ('.tfs')		# All outputs are .ptc files
madx_iterators = []				# Integers (turn) used to iterate over files

# Search from current directory
print ('\nFind MADX Twiss files\n')
for subdir, dirs, files in os.walk('MADX_Twiss'):
	# Iterate over all files
    for file in files:
		# Find files with required extension
		ext = os.path.splitext(file)[-1].lower()
		if ext in madx_extensions:
			# Print full file path
			# ~ print (os.path.join(subdir, file))		# full path to file
			fileno = int(file.split('.')[0])	# use turn number as a key
			madx_iterators.append(fileno)

madx_data = dict()

# iterate over files (turns)
print ('\nRead s, x madx_data from files\n')
last_s = 0

for turn in sorted(madx_iterators):	
	s = []
	x = []

	# Open file
	infile = 'MADX_Twiss/' + str(turn) + '.tfs'
	fin=open(infile,'r').readlines()[47:]

	# Save s, x
	for l in fin:
		if last_s == float(l.split()[2]):
			pass
		else:
			last_s =float(l.split()[2])
			s.append(float(l.split()[2]))
			x.append(float(l.split()[-2])*1E3)

	# Add to dictionary as dict[turn] = (s, x)
	madx_data[turn] = [s, x]
	last_s = 0


# Access turn 0, s column
# ptc_data[0][0]
# Access turn 25, x column, first value
# ~ print ptc_data[25][1][0]

print ('length of ptc_data = ', len(ptc_data[25][0]))
print ('max x ptc_data = ', max(ptc_data[25][1]))
print ('min x of ptc_data = ', min(ptc_data[25][1]))


#-----------------------------------------------------------------------
#------------------------------PLOTTING---------------------------------
#-----------------------------------------------------------------------

case = 'BSW_As_SBENDS'
verbose = False

print ('\n\tStart Plotting\n')

#------------------------PTC cf MADX no limits--------------------------
fig, ax1 = plt.subplots();
plt.title("PTC vs MADX Injection Closure Tune Swing");

# colormap 
colors = cm.rainbow(np.linspace(0, 1, len(ptc_iterators)))

ax1.set_xlabel("S [m]");
ax1.set_ylabel("x [m]");

c_it = int(0)
for turn in sorted(ptc_iterators):
	if verbose: print ('Plotting PTC turn ', turn)
	plt.plot(ptc_data[turn][0], ptc_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

c_it = int(0)
for turn in sorted(madx_iterators):
	if verbose: print ('Plotting MADX turn ', turn)
	plt.plot(madx_data[turn][0], madx_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                Line2D([0], [0], color=colors[-1], lw=2)]

ax1.legend(custom_lines, ['Initial Turn', 'Final Turn'])

ax1.grid(lw=0.5, ls=':');

savename = 'PTC_cf_MADX_Closed_Orbit' + case + '.png'
plt.savefig(savename, dpi = 300);

print ('\n\tPTC cf MADX no limits plot done\n')

#--------------------------PTC cf MADX zoom-----------------------------

fig, ax1 = plt.subplots();
plt.title("PTC vs MADX Injection Closure Tune Swing");

# colormap 
colors = cm.rainbow(np.linspace(0, 1, len(ptc_iterators)))

ax1.set_xlim(470.0, 510.0)

ax1.set_xlabel("S [m]");
ax1.set_ylabel("x [mm]");

test  = [1]

# ~ plt.plot(ptc_data[1][0], ptc_data[1][1], color=colors[c_it])

c_it = int(0)
for turn in sorted(ptc_iterators):
	if verbose: print ('Plotting PTC turn ', turn)
	plt.plot(ptc_data[turn][0], ptc_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

c_it = int(0)
for turn in sorted(madx_iterators):
	if verbose: print ('Plotting MADX turn ', turn)
	plt.plot(madx_data[turn][0], madx_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

ax1.vlines(BSW40, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW40, -35, 'BSW40', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW42, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW42, -35, 'BSW42', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW43, 0, 30, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW43, 25, 'BSW43', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW44, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW44, -35, 'BSW44', rotation=90, verticalalignment='center', fontsize='small')

custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                Line2D([0], [0], color=colors[-1], lw=2)]

ax1.legend(custom_lines, ['Initial Turn', 'Final Turn'])

ax1.grid(lw=0.5, ls=':');

savename = 'PTC_cf_MADX_Closed_Orbit' + case + '_zoom.png'
plt.savefig(savename, dpi = 300);

print ('\n\tPTC cf MADX zoom plot done\n')

#-----------------------------MADX zoom---------------------------------

fig, ax1 = plt.subplots();
plt.title("MADX Injection Closure Tune Swing");

# colormap 
colors = cm.rainbow(np.linspace(0, 1, len(ptc_iterators)))

ax1.set_xlim(470.0, 510.0)

ax1.set_xlabel("S [m]");
ax1.set_ylabel("x [mm]");

c_it = int(0)
for turn in sorted(madx_iterators):
	if verbose: print ('Plotting MADX turn ', turn)
	plt.plot(madx_data[turn][0], madx_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

ax1.vlines(BSW40, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW40, -35, 'BSW40', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW42, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW42, -35, 'BSW42', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW43, 0, 30, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW43, 25, 'BSW43', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW44, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW44, -35, 'BSW44', rotation=90, verticalalignment='center', fontsize='small')

custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                Line2D([0], [0], color=colors[-1], lw=2)]

ax1.legend(custom_lines, ['Initial Turn', 'Final Turn'])

ax1.grid(lw=0.5, ls=':');

plt.tight_layout()
savename = 'MADX_Closed_Orbit' + case + '_zoom.png'
plt.savefig(savename, dpi = 300);

print ('\n\tMADX zoom plot done\n')

#---------------------------MADX no limits------------------------------

fig, ax1 = plt.subplots();
plt.title("MADX Injection Closure Tune Swing");

# colormap 
colors = cm.rainbow(np.linspace(0, 1, len(ptc_iterators)))

ax1.set_xlabel("S [m]");
ax1.set_ylabel("x [mm]");

c_it = int(0)
for turn in sorted(madx_iterators):
	if verbose: print ('Plotting MADX turn ', turn)
	plt.plot(madx_data[turn][0], madx_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                Line2D([0], [0], color=colors[-1], lw=2)]

ax1.legend(custom_lines, ['Initial Turn', 'Final Turn'])

ax1.grid(lw=0.5, ls=':');

plt.tight_layout()
savename = 'MADX_Closed_Orbit' + case + '.png'
plt.savefig(savename, dpi = 300);

print ('\n\tMADX no limits plot done\n')

#--------------------------------PTC Zoom-------------------------------

fig, ax1 = plt.subplots();
plt.title("PTC Injection Closure Tune Swing");

# colormap 
colors = cm.rainbow(np.linspace(0, 1, len(ptc_iterators)))

ax1.set_xlim(470.0, 510.0)

ax1.set_xlabel("S [m]");
ax1.set_ylabel("x [mm]");

test  = [1]

# ~ plt.plot(ptc_data[1][0], ptc_data[1][1], color=colors[c_it])

c_it = int(0)
for turn in sorted(ptc_iterators):
	if verbose: print ('Plotting PTC turn ', turn)
	plt.plot(ptc_data[turn][0], ptc_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1

ax1.vlines(BSW40, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW40, -35, 'BSW40', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW42, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW42, -35, 'BSW42', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW43, 0, 30, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW43, 25, 'BSW43', rotation=90, verticalalignment='center', fontsize='small')

ax1.vlines(BSW44, -40, 0, colors='k', linestyles='solid', label='', lw=0.5, linestyle='--')
ax1.text(BSW44, -35, 'BSW44', rotation=90, verticalalignment='center', fontsize='small')

custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                Line2D([0], [0], color=colors[-1], lw=2)]

ax1.legend(custom_lines, ['Initial Turn', 'Final Turn'])

ax1.grid(lw=0.5, ls=':');

plt.tight_layout()
savename = 'PTC_Closed_Orbit' + case + '_zoom.png'
plt.savefig(savename, dpi = 300);

print ('\n\tPTC Zoom plot done\n')

#------------------------------PTC no limits----------------------------

fig, ax1 = plt.subplots();
plt.title("PTC Injection Closure Tune Swing");

# colormap 
colors = cm.rainbow(np.linspace(0, 1, len(ptc_iterators)))

ax1.set_xlabel("S [m]");
ax1.set_ylabel("x [mm]");

test  = [1]

# ~ plt.plot(ptc_data[1][0], ptc_data[1][1], color=colors[c_it])

c_it = int(0)
for turn in sorted(ptc_iterators):
	if verbose: print ('Plotting PTC turn ', turn)
	plt.plot(ptc_data[turn][0], ptc_data[turn][1], color=colors[c_it])
	# For each turn plot s,x in a new colour
	c_it += 1
        
custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                Line2D([0], [0], color=colors[-1], lw=2)]

ax1.legend(custom_lines, ['Initial Turn', 'Final Turn'])

ax1.grid(lw=0.5, ls=':');

plt.tight_layout()
savename = 'PTC_Closed_Orbit' + case + '.png'
plt.savefig(savename, dpi = 300);

print ('\n\tPTC no limits plot done\n')

#-----------------------------------------------------------------------

print ('\n\tPTC MAD-X Injection Bump Comparison Plotting Script Complete!\n')
