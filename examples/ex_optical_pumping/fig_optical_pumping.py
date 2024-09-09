import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'cm'

fig, (ax0) = plt.subplots(nrows=1, figsize=(6,6))

t, rho11_full = np.genfromtxt(r'rho11_full.dat', unpack=True)
t, rho22_full = np.genfromtxt(r'rho22_full.dat', unpack=True)
t, rho33_full = np.genfromtxt(r'rho33_full.dat', unpack=True)
t, rho44_full = np.genfromtxt(r'rho44_full.dat', unpack=True)
t, rho55_full = np.genfromtxt(r'rho55_full.dat', unpack=True)
t, rho11_rate = np.genfromtxt(r'rho11_rate.dat', unpack=True)
t, rho22_rate = np.genfromtxt(r'rho22_rate.dat', unpack=True)
t, rho33_rate = np.genfromtxt(r'rho33_rate.dat', unpack=True)
t, rho44_rate = np.genfromtxt(r'rho44_rate.dat', unpack=True)
t, rho55_rate = np.genfromtxt(r'rho55_rate.dat', unpack=True)

#  The times contained in the array t are expressed in mus. 
#  The rates were calculated for Gamma/2pi = 1 MHz. Calculating
#  Gamma t thus amounts to multiplying the content of t by 2 pi:
Gammat= (2.0*np.pi)*t

plt.plot(Gammat,rho11_full,"-",color='C{0:d}'.format(0), linewidth=1.0)
plt.plot(Gammat,rho22_full,"-",color='C{0:d}'.format(1), linewidth=1.0)
plt.plot(Gammat,rho33_full,"-",color='C{0:d}'.format(2), linewidth=1.0)
plt.plot(Gammat,rho44_full,"-",color='C{0:d}'.format(3), linewidth=1.0)
plt.plot(Gammat,rho55_full,"-",color='C{0:d}'.format(4), linewidth=1.0)

plt.plot(Gammat,rho11_rate,"--",color='C{0:d}'.format(0), linewidth=2.0)
plt.plot(Gammat,rho22_rate,"--",color='C{0:d}'.format(1), linewidth=2.0)
plt.plot(Gammat,rho33_rate,"--",color='C{0:d}'.format(2), linewidth=2.0)
plt.plot(Gammat,rho44_rate,"--",color='C{0:d}'.format(3), linewidth=2.0)
plt.plot(Gammat,rho55_rate,"--",color='C{0:d}'.format(4), linewidth=2.0)

ax0.set_xlim([-20,620])
ax0.set_ylim([-0.05,1.05])
ax0.set_xticks(np.arange(0.0, 700, 100))
ax0.set_yticks(np.arange(0.0, 1.2, 0.2))
ax0.set_xlabel(r'$\Gamma t$', fontsize = 16)
ax0.set_ylabel(r'$\rho_{ii}$', fontsize = 16)

ax0.tick_params(direction='in',bottom=True,top=True,left=True,right=True)

fig.tight_layout
plt.savefig('figure_optical_pumping.pdf', bbox_inches="tight", format='pdf')

plt.show()
