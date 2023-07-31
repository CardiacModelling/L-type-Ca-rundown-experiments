
"""
Scatter comparison of all rundown/kilosec against calcium brought in from:
Left: ICaL
Centre: Ileak
Right: ICaL + Ileak
"""
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np

import helpers

# Set font
matplotlib.rc('font', family='arial', size = 10)

temp_map = {310: 'BT', 298 : 'RT'}
color_map = {'Linear': '#1f77b4', 'Saturating': '#ff7f0e', 'Other': 'grey'}

# Index for each step
protocol = pd.read_csv('resources/protocol.csv', delimiter=',')
time = protocol.iloc[:,0]
dt = time[1] - time[0] #ms
i_start = protocol[time == 860].index.tolist()[0]
i_end = protocol[time == 1010].index.tolist()[0]

_, _, X_Ca, _, _, _, _, total = helpers.leak_proportion_calcium('BT', protocol.iloc[:, 1])
X_BT_arr = X_Ca/total
_, _, X_Ca, _, _, _, _, total = helpers.leak_proportion_calcium('BT', pd.DataFrame([-90])[0])
X_BT_hold = X_Ca/total

_, _, X_Ca, _, _, _, _, total = helpers.leak_proportion_calcium('RT', protocol.iloc[:, 1])
X_RT_arr = X_Ca/total
_, _, X_Ca, _, _, _, _, total = helpers.leak_proportion_calcium('RT', pd.DataFrame([-90])[0])
X_RT_hold = X_Ca/total

leak_prop = {'BT': [X_BT_arr, X_BT_hold], 'RT': [X_RT_arr, X_RT_hold]}

# load rrate, shape, cell name
rrate_data = pd.read_csv('output/r_rate_database.csv')

Ca_frac = []
Ca_eff = []
temp_arr = []
thold_arr = []
inaca_arr = []
Ca_ntot = []

fig = plt.figure(figsize=(6.6, 3))
ax_ical = fig.add_subplot(131)
ax_ileak = fig.add_subplot(132)
ax_itot = fig.add_subplot(133)

for i in range(len(rrate_data)):
    cell = rrate_data['Cell ID'].iloc[i]
    r_rate = rrate_data['Run rate'].iloc[i]
    shape = rrate_data['shape'].iloc[i]

    temp = temp_map[rrate_data['Temperature'].iloc[i]]
    hold = rrate_data['thold'].iloc[i]

    temp_arr.append(rrate_data['Temperature'].iloc[i])
    inaca_arr.append(rrate_data['INaCa'].iloc[i])
    # print(rrate_data['Temperature'].iloc[i])
    # continue

    # load the gleak, Eleak, and Cap 
    pathtoprop = f'output/{temp}_{hold}/prop_{cell}.csv'
    prop_data = pd.read_csv(pathtoprop)

    gleak = prop_data['gleak (nS)']
    Eleak = prop_data['Eleak (mV)']
    cap = prop_data['cap (pF)']

    # load the ical across all sweeps for the cell
    pathtocell = f'output/{temp}_{hold}/{cell}.csv'
    ical_all = pd.read_csv(pathtocell)

    n_sweeps = len(ical_all.columns)

    # Calculate NCa from ical
    Cai = - dt *ical_all.iloc[i_start: i_end, :].sum(axis = 0)/ (2 * 96500) # fmol
    Cai_ical = Cai.sum()

    # Calculate IleakCa
    ## during sweeps
    ileak_swe = pd.DataFrame()
    for i in range(n_sweeps):
        ileak_swe[i] = gleak.iloc[i] * (protocol.iloc[:, 1] - Eleak.iloc[i]) * leak_prop[temp][0]
    Ca_leak = - dt *ileak_swe.iloc[i_start: i_end, :].sum(axis = 0)/ (2 * 96500) # fmol
    Ca_leak = Ca_leak.sum()

    ## between sweeps
    t_sweep = pd.read_csv(f'resources/{temp}_{hold}_sweep_time.csv')
    for i in range(n_sweeps - 1):
        if i == 0:
            I = gleak.iloc[i] * (-90 - Eleak.iloc[i]) * leak_prop[temp][1]
            ca = - 1000 * (t_sweep.iloc[i+1] - 0)*I[0]/(2 * 96500) # fmol 
            Ca_leak += ca[0]
        elif np.isnan(gleak.iloc[i]):
            I = gleak.iloc[i-1] * (-90 - Eleak.iloc[i-1]) * leak_prop[temp][1]
            ca = - 1000 * (t_sweep.iloc[i+1] - t_sweep.iloc[i])*I[0]/(2 * 96500) # fmol 
            Ca_leak += ca[0]
        else:
            I = gleak.iloc[i] * (-90 - Eleak.iloc[i]) * leak_prop[temp][1]
            ca = - 1000 * (t_sweep.iloc[i+1] - t_sweep.iloc[i])*I[0]/(2 * 96500) # fmol 
            Ca_leak += ca[0]

    # Normalise NCa with Cap
    cap_area = cap.median()**1.5
    Ca_ical_norm = Cai_ical/cap_area
    Ca_leak_norm = Ca_leak/cap_area
    Ca_tot_norm = Ca_ical_norm + Ca_leak_norm

    Ca_frac.append(Ca_leak/(Ca_leak + Cai_ical))
    Ca_eff.append(Ca_tot_norm)
    thold_arr.append(hold)
    Ca_ntot.append(Ca_leak + Cai_ical)

    # plot
    ax_ical.scatter(Ca_ical_norm, r_rate, edgecolors = color_map[shape], lw = 2, facecolors = 'none' )
    ax_ileak.scatter(Ca_leak_norm, r_rate, edgecolors = color_map[shape], lw = 2, facecolors = 'none' )
    ax_itot.scatter(Ca_tot_norm, r_rate, edgecolors = color_map[shape], lw = 2, facecolors = 'none' )

# Ca_frac = pd.DataFrame(Ca_frac, columns=['leak moles frac', 'Ca_eff'])
df = {'leak moles frac': Ca_frac, 'Ca_eff': Ca_eff, 'Ca_tot': Ca_ntot, 'Temperature': temp_arr, \
      'thold': thold_arr, 'INaCa': inaca_arr}
df = pd.DataFrame(df)
df.to_csv(f'resources/ca_leak_frac.csv')

ax_ileak.scatter([], [], edgecolors = color_map['Linear'], lw = 2, facecolors = 'none', label = 'linear')
ax_ileak.scatter([], [], edgecolors = color_map['Saturating'], lw = 2, facecolors = 'none', label = 'saturating')
ax_ileak.scatter([], [], edgecolors = color_map['Other'], lw = 2, facecolors = 'none', label = 'other' )

ax_ical.set_ylabel('$R_{rate}$ (1/min)')
ax_ical.set_xlabel('$N_{Ca}$/$C_{m}^{{3}/{2}}$ (fmol/$\mathrm{pF}^{3/2}$)')
# ax_ical.set_ylim(0, 3)
ax_ical.set_xlim(0, 0.2)


ax_ileak.set_xlabel('$N_{Ca}$/$C_{m}^{{3}/{2}}$ (fmol/$\mathrm{pF}^{3/2}$)')
# ax_ileak.set_ylim(0, 3)
ax_ileak.set_xlim(0, 0.11)
ax_ileak.legend(ncol =3, bbox_to_anchor = (0.2, 1.02, 0.6, 0.15), loc = 'center', bbox_transform = ax_ileak.transAxes)

ax_itot.set_xlabel('$N_{Ca}$/$C_{m}^{{3}/{2}}$ (fmol/$\mathrm{pF}^{3/2}$)')
# ax_itot.set_ylim(0, 3)
ax_itot.set_xlim(0, 0.28)


plt.subplots_adjust(bottom=0.17, wspace = 0.3)
plt.savefig('figures/figure5.pdf')
plt.close()