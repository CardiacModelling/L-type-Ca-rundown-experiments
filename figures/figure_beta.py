import matplotlib
import matplotlib.pyplot as plt

# Set font
matplotlib.rc('font', family='arial', size = 10)

# defaults: INaCa_off = 0, INaCa_on = 1; temp_310 = 0, temp_298 = 1
# equation: r_rate = b_0 + b_1 * (log(thold)) + b_2 * (temp_298) + b_3 * (INaCA_On)
b_0, b_1, b_2, b_3 = 0.158795, -0.016373, -0.062213, -0.021030
ci_l_0, ci_u_0 =  0.11984411,  0.197745292
ci_l_1, ci_u_1 = -0.02930387, -0.003442791
ci_l_2, ci_u_2 = -0.07628548, -0.048139628
ci_l_3, ci_u_3 = -0.03490557, -0.007153578

# Plot Confidence interval and beta values
fig = plt.figure(figsize=(6.6, 3))

plt.scatter(0, -b_0, color = 'tab:blue')
plt.hlines(- ci_l_0, xmin = -0.2, xmax = 0.2, color = 'grey')
plt.hlines(- ci_u_0, xmin = -0.2, xmax = 0.2, color = 'grey')
plt.vlines(0, ymin = - ci_l_0, ymax = - ci_u_0, color = 'grey')

plt.scatter(1, b_1, color = 'tab:blue')
plt.hlines(ci_l_1, xmin = 0.8, xmax = 1.2, color = 'grey')
plt.hlines(ci_u_1, xmin = 0.8, xmax = 1.2, color = 'grey')
plt.vlines(1, ymin = ci_l_1, ymax = ci_u_1, color = 'grey')

plt.scatter(2, b_2, color = 'tab:blue')
plt.hlines(ci_l_2, xmin = 1.8, xmax = 2.2, color = 'grey')
plt.hlines(ci_u_2, xmin = 1.8, xmax = 2.2, color = 'grey')
plt.vlines(2, ymin = ci_l_2, ymax = ci_u_2, color = 'grey')

plt.scatter(3, b_3, color = 'tab:blue')
plt.hlines(ci_l_3, xmin = 2.8, xmax = 3.2, color = 'grey')
plt.hlines(ci_u_3, xmin = 2.8, xmax = 3.2, color = 'grey')
plt.vlines(3, ymin = ci_l_3, ymax = ci_u_3, color = 'grey')

plt.xticks([0, 1, 2, 3], [r'-$\beta_0$', r'$\beta_1$', r'$\beta_2$', r'$\beta_3$'])
plt.ylabel('Coefficient')

plt.tight_layout()
plt.savefig('figures/figure_beta.pdf', facecolor='w', transparent=False)
plt.close()