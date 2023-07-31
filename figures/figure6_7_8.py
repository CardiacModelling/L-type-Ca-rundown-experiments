"""
Plot rundown rate split by experiemntal condition
"""
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

matplotlib.rc('font', family='arial', size = 10)

fig = [plt.figure(figsize=(6.6, 3.2)) for i in range(3)]
sub = [fig[i].add_subplot(111) for i in range(3)]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
sc_size = 0.06


# load rrate, shape, cell name
data = pd.read_csv('output/r_rate_database.csv')

# temp = BT, thold = 10, INaCa = On
rate = data[(data['Temperature'] == 310) & (data['thold'] == 10) & (data['INaCa'] == 'On')]['Run rate']
sub[0].boxplot(rate, positions = [0.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(0.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[0], label = '$t_{hold}$: 10s', lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [0.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
x = np.random.normal(0.5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[0], label = '310 K', lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [1], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
x = np.random.normal(1, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

# temp = BT, thold = 10, INaCa = Off
rate = data[(data['Temperature'] == 310) & (data['thold'] == 10) & (data['INaCa'] == 'Off')]['Run rate']
sub[0].boxplot(rate, positions = [2.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(2.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [2], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(2, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [0.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(0.5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[0], label = 'INaCa Blocked', lw = 2, facecolors = 'none')

# temp = BT, thold = 20, INaCa = On
rate = data[(data['Temperature'] == 310) & (data['thold'] == 20) & (data['INaCa'] == 'On')]['Run rate']
sub[0].boxplot(rate, positions = [1], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(1, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[1], label = '$t_{hold}$: 20s', lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [3.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(3.5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [2.5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(2.5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none', label = 'No Block')

# temp = BT, thold = 20, INaCa = Off
rate = data[(data['Temperature'] == 310) & (data['thold'] == 20) & (data['INaCa'] == 'Off')]['Run rate']
sub[0].boxplot(rate, positions = [3], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(3, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [2], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(2, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

# temp = BT, thold = 40, INaCa = On
rate = data[(data['Temperature'] == 310) & (data['thold'] == 40) & (data['INaCa'] == 'On')]['Run rate']
sub[0].boxplot(rate, positions = [1.5], medianprops = \
    dict(color = colors[2]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(1.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[2], label = '$t_{hold}$: 40s', lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [6.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(6.5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [4], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(4, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

# temp = BT, thold = 40, INaCa = Off
rate = data[(data['Temperature'] == 310) & (data['thold'] == 40) & (data['INaCa'] == 'Off')]['Run rate']
sub[0].boxplot(rate, positions = [3.5], medianprops = \
    dict(color = colors[2]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(3.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[2], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [8], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(8, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [3.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(3.5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

# temp = RT, thold = 10, INaCa = On
rate = data[(data['Temperature'] == 298) & (data['thold'] == 10) & (data['INaCa'] == 'On')]['Run rate']
sub[0].boxplot(rate, positions = [4.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(4.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [1], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(1, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[1], label = '298 K', lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [5.5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(5.5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

# temp = RT, thold = 10, INaCa = Off
rate = data[(data['Temperature'] == 298) & (data['thold'] == 10) & (data['INaCa'] == 'Off')]['Run rate']
sub[0].boxplot(rate, positions = [6.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(6.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [2.5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(2.5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')


# temp = RT, thold = 20, INaCa = On
rate = data[(data['Temperature'] == 298) & (data['thold'] == 20) & (data['INaCa'] == 'On')]['Run rate']
sub[0].boxplot(rate, positions = [5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [4], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(4, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [7], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(7, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

# temp = RT, thold = 10, INaCa = Off
rate = data[(data['Temperature'] == 298) & (data['thold'] == 20) & (data['INaCa'] == 'Off')]['Run rate']
sub[0].boxplot(rate, positions = [7], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(7, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [5.5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(5.5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [6.5], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(6.5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')


# temp = RT, thold = 10, INaCa = On
rate = data[(data['Temperature'] == 298) & (data['thold'] == 40) & (data['INaCa'] == 'On')]['Run rate']
sub[0].boxplot(rate, positions = [5.5], medianprops = \
    dict(color = colors[2]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(5.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[2], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [7], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(7, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [8.5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(8.5, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

# temp = RT, thold = 40, INaCa = On
rate = data[(data['Temperature'] == 298) & (data['thold'] == 40) & (data['INaCa'] == 'Off')]['Run rate']
sub[0].boxplot(rate, positions = [7.5], medianprops = \
    dict(color = colors[2]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(7.5, 0.08, size = len(y))
sub[0].scatter(x, y, edgecolors = colors[2], lw = 2, facecolors = 'none')

sub[1].boxplot(rate, positions = [8.5], medianprops = \
    dict(color = colors[1]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(8.5, 0.08, size = len(y))
sub[1].scatter(x, y, edgecolors = colors[1], lw = 2, facecolors = 'none')

sub[2].boxplot(rate, positions = [8], medianprops = \
    dict(color = colors[0]), widths = 0.4, showfliers = False, whis = [0, 100])
y = rate.dropna()
x = np.random.normal(8, 0.08, size = len(y))
sub[2].scatter(x, y, edgecolors = colors[0], lw = 2, facecolors = 'none')

sub[0].set_xticks([1, 3, 5, 7])
sub[1].set_xticks([0.75, 2.25, 3.75, 5.25, 6.75, 8.25])
sub[2].set_xticks([0.75, 2.25, 3.75, 5.25, 6.75, 8.25])

index_hold = ['310K\nNo Block', '310K\nINaCa Blocked', '298K\nNo Block', '298K\nINaCa Blocked']
sub[0].set_xticklabels(index_hold)

index_temp = ['$t_{hold}$: 10s\nNo Block', '$t_{hold}$: 10s\nINaCa Blocked',\
    '$t_{hold}$: 20s\nNo Block', '$t_{hold}$: 20s\nINaCa Blocked',\
    '$t_{hold}$: 40s\nNo Block', '$t_{hold}$: 40s\nINaCa Blocked']
sub[1].set_xticklabels(index_temp)

index_inaca = ['$t_{hold}$: 10s\n310K', '$t_{hold}$: 20s\n310K',\
    '$t_{hold}$: 40s\n310K', '$t_{hold}$: 10s\n298K',\
    '$t_{hold}$: 20s\n298K', '$t_{hold}$: 40s\n298K']
sub[2].set_xticklabels(index_inaca)


for el in sub:
    el.set_ylabel('$R_{rate}$ (1/min)')
    #el.set_ylim(-0.003, 0.004)
    el.legend(loc = 'lower left')

fig[0].savefig('figures/figure6.pdf')
fig[1].savefig('figures/figure8.pdf')
fig[2].savefig('figures/figure7.pdf')

plt.tight_layout()
plt.close()