import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

# Set font
matplotlib.rc('font', family='arial', size = 10)

# create data
protocol_fig = pd.read_csv('resources/protocol_fig.csv', delimiter=',')
time = protocol_fig.iloc[:,0].to_numpy()/1000
volt = protocol_fig.iloc[:,1].to_numpy()

x = np.concatenate((time, time + 10))
y = np.concatenate((volt, volt))

# Create two subplots with shared y-axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (6.6, 3))

# Set the limits of the x-axis on each subplot
ax1_x1 = -0.2
ax1_x2 = 1.2

ax2_x1 = 9.8
ax2_x2 = 11.2

ax1.set_xlim(ax1_x1, ax1_x2)
ax2.set_xlim(ax2_x1, ax2_x2)

y_min = min(y) - 5
y_max = max(y) + 5

ax1.set_ylim(y_min, y_max)
ax2.set_ylim(y_min, y_max)

ax1.text(ax1_x2 - 0.01, y_min - 5, '/', size = 15)
ax2.text(ax2_x1 - 0.02, y_min - 5, '/', size = 15)

ax1.text(ax1_x2 - 0.01, y_max - 5, '/', size = 15)
ax2.text(ax2_x1 - 0.02, y_max - 5, '/', size = 15)

ax1.text(ax1_x2 - 0.01, -90 - 5, '/', size = 15)
ax2.text(ax2_x1 - 0.02, -90 - 5, '/', size = 15)

# Plot the data on each subplot
ax1.plot(x[x<10], y[x<10], 'tab:blue')
ax2.plot(x[x>=10], y[x>=10], 'tab:blue')

# plot fillers
ax1.plot([ax1_x1 + 0.05, 0], [-90, -90], 'tab:blue')
ax1.plot([max(time), ax1_x2], [-90, -90], 'tab:blue')

ax2.plot([ax2_x1, 10], [-90, -90], 'tab:blue')
ax2.plot([10 + max(time), ax2_x2 - 0.05], [-90, -90], 'tab:blue')

ax2.set_yticks([])

# Hide the spines (the lines connecting the axis tick marks) on the right of ax1 and on the left of ax2
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)

# label the times
ax1.text(ax1_x2, -85, r'$\mathrm{t}$=')
ax2.text(ax2_x1, -85, r'$\mathrm{t}_{hold}$')

ax1.set_ylabel('Voltage (mV)')
ax1.set_xlabel('Time (s)', horizontalalignment='right', x=1.15)

plt.tight_layout(w_pad = 0.01)

plt.savefig('figures/figure_prot.pdf')
plt.close()