# Creates a figure showing post-processing of the raw current
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
# import numpy as np
import os

import sys
sys.path.append('../')
import helpers # for cap ciltering

# Set font
matplotlib.rc('font', family='arial', size = 10)

# load protocol
protocol = pd.read_csv('../resources/protocol.csv', delimiter=',')
time = protocol.iloc[:,0]
ind_step_start = protocol[time == 860].index.tolist()[0]
ind_step_end = protocol[time == 1010].index.tolist()[0]

# load all current
directoryname = 'Cav1.2_Run_Down_BT_10s_20210129_12.07.03'
pathtodatadirectory = 'raw_data/'
files = os.listdir(f'../raw_data/{directoryname}')
cell = 'C19'
for f in files:
    if f[-8:-5] == cell:
        filename = f
        path = f'../raw_data/{directoryname}/{filename}'
        break
cell_data = pd.read_csv(path, sep =';', usecols=range(1,58 + 2))

all_sweeps = cell_data.iloc[2:, 1:].astype(float) * pow(10, 12) # Convert A to pA 

all_sweeps = helpers.cap_filter(ind_step_start, ind_step_end, 0.1, all_sweeps) #Cap filter
all_sweeps = all_sweeps.iloc[:, :34 + 1]

# load the gleak, Eleak, etc.
prop_data = pd.read_csv(f'../output/BT_10/prop_{cell}.csv')
g_1, E_1 = prop_data.iloc[0]['gleak (nS)'], prop_data.iloc[0]['Eleak (mV)']
g_n, E_n = prop_data.iloc[-1]['gleak (nS)'], prop_data.iloc[-1]['Eleak (mV)'] 

leak_1 = g_1*(protocol.iloc[:, 1] - E_1)
leak_2 = g_n*(protocol.iloc[:, 1] - E_n)

# leak subtracted current
sweep1_sub = all_sweeps.iloc[:, 0].reset_index(drop=True).sub(leak_1.reset_index(drop=True), axis = 0) 
sweepn_sub = all_sweeps.iloc[:, -1].reset_index(drop=True).sub(leak_2.reset_index(drop=True), axis = 0) 

# load figure
fig = plt.figure(figsize=(6.6, 8))
gs = gridspec.GridSpec(ncols =1, nrows = 5)
sub1 = fig.add_subplot(gs[0, 0])
sub2 = fig.add_subplot(gs[1, 0])
sub3 = fig.add_subplot(gs[2, 0])
sub4 = fig.add_subplot(gs[3, 0])
sub5 = fig.add_subplot(gs[4, :])

props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
sub1.text(0.015, 0.12, 'A', transform = sub1.transAxes, weight = 'bold', bbox=props)
sub2.text(0.015, 0.12, 'B', transform = sub2.transAxes, weight = 'bold', bbox=props)
sub3.text(0.015, 0.12, 'C', transform = sub3.transAxes, weight = 'bold', bbox=props)
sub4.text(0.015, 0.12, 'D', transform = sub4.transAxes, weight = 'bold', bbox=props)
sub5.text(0.015, 0.12, 'E = C-D', transform = sub5.transAxes, weight = 'bold', bbox=props)

sub1.axhline(0, color='grey', ls= 'dashed', lw = 0.8)
sub2.axhline(0, color='grey', ls= 'dashed', lw = 0.8)
sub3.axhline(0, color='grey', ls= 'dashed', lw = 0.8)
sub4.axhline(0, color='grey', ls= 'dashed', lw = 0.8)
sub5.axhline(0, color='grey', ls= 'dashed', lw = 0.8)

# plot inset protocol
protocol_fig = pd.read_csv('../resources/protocol_fig.csv', delimiter=',')
ins_ax = inset_axes(sub3, width = "30%", height = '52%', loc = 'center', borderpad =3.2)
ins_ax.plot(time, protocol_fig.iloc[:,1], color ='grey')

colors = ['#1f77b4', '#ff7f0e']

sub1.plot(time, all_sweeps.iloc[:, 0], color = colors[0])
sub2.plot(time, all_sweeps.iloc[:, -1], color = colors[1])
sub1.plot(time, leak_1, color = 'grey', label = \
          f'gleak = {round(g_1, 1)}nS, Eleak = {round(E_1, 1)}mV')
sub2.plot(time, leak_2, color = 'grey', label = \
          f'gleak = {round(g_n, 1)}nS, Eleak = {round(E_n, 1)}mV')
sub3.plot(time, sweep1_sub, color = colors[0])
sub4.plot(time, sweepn_sub, color = colors[1])
sub5.plot(time, sweep1_sub - sweepn_sub.values, color = colors[0])


sub1.set_ylabel('Raw current (pA)\n sweep 1')
sub3.set_ylabel('Leak-subtracted\ncurrent (pA)\nsweep 1')
sub5.set_ylabel('Leak-drug-\nsubtracted\ncurrent (pA)\nsweep 1')
sub2.set_ylabel('Raw current (pA)\nsweep n2')
sub4.set_ylabel('Leak-subtracted\ncurrent (pA)\nsweep n2')
sub5.set_xlabel('Time from the beginning of the first sweep (ms)')

sub2.legend(loc = 'upper center')
sub1.legend(loc = 'lower center')

yticks = matplotlib.ticker.MaxNLocator(4)
sub1.yaxis.set_major_locator(yticks)
sub3.yaxis.set_major_locator(yticks)
sub4.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
sub5.yaxis.set_major_locator(yticks)


plt.tight_layout()
plt.savefig('figure3.pdf')
plt.close()