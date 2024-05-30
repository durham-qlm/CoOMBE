import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'cm'

fig, (ax0) = plt.subplots(nrows=1,figsize=(6.4,4.0))

Delta, chi_re2, chi_im2, refindx2, alpha2 = \
                              np.genfromtxt(r'chi_weakprb_2', unpack=True)
Delta, chi_re3, chi_im3, refindx3, alpha3 = \
                              np.genfromtxt(r'chi_weakprb_3', unpack=True)
Delta, chi_renw, chi_imnw, refindxnw, alphanw = \
                              np.genfromtxt(r'chi_noweakprb', unpack=True)
plt.plot(Delta,alpha2)
plt.plot(Delta,alpha3,'--')
plt.plot(Delta,alphanw,':')


ax0.set_xlim([-2000.,2000.])
ax0.set_xticks(np.arange(-2000., 3000., 1000.))
ax0.set_yscale('log')
#ax0.set_ylim([0.0024,0.0064])
#ax0.set_yticks(np.arange(0.003, 0.007, 0.001))
ax0.set_xlabel(r'$\Delta/(2\pi)$ $\mathrm{(MHz)}$', fontsize=16)
ax0.set_ylabel(r'$\alpha$ $\mathrm{(m}^{-1}\mathrm{)}$', fontsize = 16)


ax0.tick_params(direction='in',bottom=True,top=True,left=True,right=True)

fig.tight_layout
plt.savefig('figure_Rb85D1.pdf', bbox_inches="tight", format='pdf')

plt.show()
