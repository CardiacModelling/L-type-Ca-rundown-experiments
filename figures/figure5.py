# Plot rundown for all selected cells
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import extensions
import helpers

# Index for each step
protocol = pd.read_csv('resources/protocol.csv', delimiter=',')
time = protocol.iloc[:,0]
i_start = protocol[time == 860].index.tolist()[0]
i_end = protocol[time == 1010].index.tolist()[0]


temperature = ['BT', 'RT']
holding_duration = [10, 20, 40]

# plot

# Create a 4x3 plot
fig = plt.figure(figsize=(6.6, 8))
gs = fig.add_gridspec(4, 3)

sub = [[], [], [], []]
for i in range(4):
    for j in range(3):
        sub[i].append(fig.add_subplot(gs[i, j]))
        sub[i][j].set_ylim(-0.2, 1)
        if j == 0:
            if i == 0:
                sub[i][j].set_ylabel('310 K, INaCa On\nRundown')
                sub[i][j].set_title('$t_{hold}$: 10s')
            elif i == 1:
                sub[i][j].set_ylabel('310 K, INaCa Off\nRundown')
            elif i == 2:
                sub[i][j].set_ylabel('298 K, INaCa On\nRundown')
            else:
                sub[i][j].set_ylabel('298 K, INaCa Off\nRundown')
        if i == 0:
            if j == 1:
                sub[i][j].set_title('$t_{hold}$: 20s')
            elif j == 2:
                sub[i][j].set_title('$t_{hold}$: 40s')
        if i == 3:
            sub[i][j].set_xlabel('Time (s)')

        if j == 1 or j == 2:
            ticks = [1, 0.75, 0.5, 0.25, 0] 
            sub[i][j].set_yticks(ticks = ticks, labels = [])

        if i == 0 or i == 1 or i == 2:
            ticks = [0, 100, 200, 300]
            sub[i][j].set_xticks(ticks = ticks, labels = [])

# store plots by the relevant experimental conditions
plots_dicts = {
    '10' : {
        'BT' : [sub[0][0], sub[1][0]],
        'RT' : [sub[2][0], sub[3][0]]
    },
    '20' : {
        'BT' : [sub[0][1], sub[1][1]],
        'RT' : [sub[2][1], sub[3][1]]
    },
    '40' : {
        'BT' : [sub[0][2], sub[1][2]],
        'RT' : [sub[2][2], sub[3][2]]
    }
}

for temp in temperature:
    for hold_dur in holding_duration:
        # load the data
        cells = extensions.selected_cells(temp, hold_dur)
        sweep_time = pd.read_csv(f'resources/{temp}_{hold_dur}_sweep_time.csv')
        
        path = f'output/{temp}_{hold_dur}/'

        count_on = 0
        count_off = 0
        plots = plots_dicts[f'{hold_dur}'][temp]

        for cell in cells:
            ptocell = path + f'{cell}.csv'
            ical_all = pd.read_csv(ptocell)

            # Calculate rundown array 
            ical_peak = ical_all.iloc[i_start: i_end, :].min(axis = 0)
            rundown = 1-1*ical_peak/ical_peak.iloc[0]
            rundown = rundown.tolist()

            # calculate time array 
            i_min = ical_all.iloc[i_start: i_end, :].idxmin(axis = 0)
            t_stamp = [] # in seconds
            for i in range(len(i_min)):
                if np.isnan(i_min[i]):
                    t_stamp.append(np.nan)
                else:
                    # add time from sweep to the sweep start time
                    t = protocol.iloc[int(i_min[i]), 0]/1000 + \
                        sweep_time.iloc[i].to_numpy()[0]  # seconds

                    t_stamp.append(t) # in seconds 

            t_stamp = [x for x in t_stamp if x == x] # cleans out float NaNs
            rundown = [x for x in rundown if x == x] # cleans out float NaNs
           
            # count and plot in the relevant subplot
            inaca_stat = helpers.inaca_status(cell)
            if inaca_stat == 'On':
                count_on += 1
                plots[0].plot(t_stamp, rundown, marker = 'o', color = 'grey', markersize = 0.5, lw = 0.5)
            if inaca_stat == 'Off':
                count_off += 1
                plots[1].plot(t_stamp, rundown, marker = 'o', color = 'grey', markersize = 0.5, lw = 0.5)
        
        plots[0].text(0.05, 0.9, str(count_on), transform = plots[0].transAxes)
        plots[1].text(0.05, 0.9, str(count_off), transform = plots[1].transAxes)


plt.tight_layout()
fig.savefig('figures/figure5.pdf', facecolor='w', transparent=False)
plt.close()
