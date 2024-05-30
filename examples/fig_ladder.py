import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'cm'

fig, (ax0) = plt.subplots(nrows=1,figsize=(6.4,4.0))

Delta, rho21_re, rho21_im = np.genfromtxt(r'ladder_rho21.dat', unpack=True)

plt.plot(Delta,-rho21_im)

ax0.set_xlim([-30,30])
ax0.set_xticks(np.arange(-30, 40, 10.0))
ax0.set_ylim([0.0024,0.0066])
ax0.set_yticks(np.arange(0.003, 0.007, 0.001))
ax0.set_xlabel(r'$\Delta/(2\pi)$ $\mathrm{(MHz)}$', fontsize=16)
ax0.set_ylabel(r'$-\mathrm{Im}\,\rho_{21}$', fontsize = 16)


ax0.tick_params(direction='in',bottom=True,top=True,left=True,right=True)

fig.tight_layout
plt.savefig('figure_ladder.pdf', bbox_inches="tight", format='pdf')

plt.show()
