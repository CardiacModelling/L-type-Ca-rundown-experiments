# rundown illustration on one cell BT 10, C19
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import numpy as np

# Set font
matplotlib.rc('font', family='arial', size = 10)

protocol = pd.read_csv('resources/protocol.csv', delimiter=',')
ind_step_start = protocol[protocol.iloc[:,0] == 860].\
    index.tolist()[0]
ind_step_end = protocol[protocol.iloc[:,0] == 1010].\
    index.tolist()[0]
time = protocol.iloc[:,0].iloc[ind_step_start: ind_step_end] 
time_x = time - time.iloc[0]

# Load sweep times
sweep_time = pd.read_csv('resources/BT_10_sweep_time.csv') # in seconds


data = pd.read_csv('output/BT_10/C19.csv')
fig = plt.figure(figsize=(6.6, 3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax3 = ax2.twinx()

n_col = len(data.columns)

colors = [matplotlib.cm.viridis(x) for x in np.linspace(0, 1, n_col)]

for i in range(n_col):
    ax1.plot(time_x, data.iloc[ind_step_start:ind_step_end, i], color = colors[i])
    ical_peak = data.iloc[ind_step_start:ind_step_end, i].min(axis=0)
    ind_peak = data.iloc[ind_step_start:ind_step_end, i].idxmin(axis=0) - ind_step_start
    try: # needed to exclude nan sweeps
        # calculate the time stamp from the beginning of the experiment
        t_sweep = (time_x.iloc[ind_peak] + 860)/1000 # time of peak
        t_tot = sweep_time.iloc[i][0] # time before sweep
        t = t_sweep + t_tot  # time stamp (in seconds)
        ax2.scatter(t, ical_peak, color = colors[i])
        ax3.scatter(t, 1 - ical_peak/data.iloc[ind_step_start:ind_step_end, 0].min(axis=0), colors[i])
    except:
        pass

# inset protocol
ins_ax = inset_axes(ax1, width = "30%", height = '30%', loc = 'lower right', borderpad =1.8)
ins_ax.plot(protocol.iloc[ind_step_start - 200: ind_step_end + 500, 0] - \
            protocol.iloc[ind_step_start: ind_step_end + 700, 0].iloc[0], \
    protocol.iloc[ind_step_start - 200: ind_step_end + 500,1], color ='grey')
ins_ax.set_title('X 33 times')

ax1.set_xlabel('Time from the beginning of the\n step to 0 mV for each sweep (ms)')
ax2.set_xlabel('Time (s)')
ax1.set_ylabel('Current (pA)')
ax2.set_ylabel('Peak Current (pA)')
norm_f = data.iloc[ind_step_start:ind_step_end, 0].min(axis=0)
ymin, ymax = ax1.get_ylim()
ymin_2, ymax_2 = ax2.get_ylim()
ymin_3, ymax_3 = ax3.get_ylim()
ax2.set_ylim(ymin, ymax)
ax3.set_ylim(1 - ymin/norm_f, 1-ymax/norm_f)
ax3.set_ylabel('Rundown')
plt.tight_layout()
plt.savefig('figures/figure1.pdf')
