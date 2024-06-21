import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'cm'

fig, (ax0) = plt.subplots(nrows=1,figsize=(6.4,4.0))

t, rho11 = np.genfromtxt(r'timedep_rho11.dat', unpack=True)
t, rho12_re, rho12_im = np.genfromtxt(r'timedep_rho12.dat', unpack=True)
t = t*1000 # Conversion from mus to ns.

plt.plot(t,rho11)
plt.plot(t,rho12_im,"--")

ax0.set_xlim([-2,2])
ax0.set_xticks(np.arange(-2.0, 3.0, 1.0))
ax0.set_yticks(np.arange(-0.5, 1.5, 0.5))
ax0.set_xlabel(r'$t$ $\mathrm{(ns)}$', fontsize=16)
ax0.set_ylabel(r'$\rho_{11},\;\mathrm{Im}\,\rho_{12}$', fontsize = 16)

ax0.tick_params(direction='in',bottom=True,top=True,left=True,right=True)

fig.tight_layout
plt.savefig('figure_timedep.pdf', bbox_inches="tight", format='pdf')

plt.show()
