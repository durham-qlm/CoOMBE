import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'cm'

cmap_steps = 256 # number of steps on colour map

N = 150 # number of points for plotting/interpolation

fig, (ax0,ax1) = plt.subplots(nrows=2,sharex=True,figsize=(6.4,6.0))

t, z, probe_re, probe_im, coupl_re, coupl_im = np.genfromtxt(r'outamplitudes.dat', unpack=True)
t = t*1000 # Conversion from mus to ns.
ti = np.linspace(-1.0, t.max(), N)
zi = np.linspace(z.min(), z.max(), N)
#  Conversion from field strenghts to intensity:
factor = 10.e-6 / probe_re[0]**2  # 10 muW/cm2.
probe = factor*(probe_re**2 + probe_im**2)
coupl = factor*(coupl_re**2 + coupl_im**2)

################################################################

probe_interp = scipy.interpolate.griddata((t, z), probe, (ti[None,:], zi[:,None]), method='linear')
cmap_range_probe = np.linspace(0., probe_interp.max(), cmap_steps)
cf0 = ax0.contourf(ti,zi,probe_interp,cmap_range_probe, cmap=plt.cm.gist_heat_r)
cbar0 = fig.colorbar(cf0,ax=ax0,ticks=[0,6e-5])
cbar0.ax.set_yticklabels(['0.','60.'])
cbar0.set_label('$\mathrm{Intensity}$ $\mathrm{(}\mu\mathrm{W}$ $\mathrm{cm}^{-2}\mathrm{)}$', fontsize = 16)
ax0.set_xlim([-1,2])
ax0.set_xticks(np.arange(-1.0, 2.5, 0.5))
ax0.set_yticks(np.arange(0, 20, 4))
ax0.set_ylabel(r'$z$ $\mathrm{(}\mu\mathrm{m)}$', fontsize = 16)
ax0.tick_params(direction='in',bottom=True,top=True,left=True,right=True)
ax0.text(-0.9,13.8,'$\mathrm{(a)}$',fontsize = 18)

################################################################

coupl_interp = scipy.interpolate.griddata((t, z), coupl, (ti[None,:], zi[:,None]), method='linear')
cmap_range_coupl = np.linspace(0., coupl_interp.max(), cmap_steps)
cf1 = ax1.contourf(ti,zi,coupl_interp,cmap_range_coupl, cmap=plt.cm.Greys)
cbar1 = fig.colorbar(cf1,ax=ax1,ticks=[0,2.5e3])
cbar1.ax.set_yticklabels(['0.0','2.5'])
cbar1.set_label('$\mathrm{Intensity}$ $\mathrm{(}\mathrm{kW}$ $\mathrm{cm}^{-2}\mathrm{)}$', fontsize = 16)
ax1.set_xlim([-1,2])
ax1.set_xticks(np.arange(-1.0, 2.5, 0.5))
ax1.set_yticks(np.arange(0, 20, 4))
ax1.set_ylabel(r'$z$ $\mathrm{(}\mu\mathrm{m)}$', fontsize = 16)
ax1.tick_params(direction='in',bottom=True,top=True,left=True,right=True)
ax1.text(-0.9,13.8,'$\mathrm{(b)}$',fontsize = 18)
ax1.set_xlabel(r'$t - z/c$ $\mathrm{(ns)}$', fontsize=16)

################################################################

fig.tight_layout()

plt.subplots_adjust(hspace=0.115)

plt.savefig('figure_prop.pdf', format='pdf')

plt.show()
