# ::setlocal makeprg=cd\ script\ &&\ python\ plots.py
import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 13
MEDIUM_SIZE = 18
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize


def plot_Vs(l_per_M, save=['yes','no']):

    r_per_M = np.arange(2,40,0.1)
    Veff = l_per_M**2 * r_per_M**-2 / 2 - r_per_M**-1 - l_per_M**2 * r_per_M**-3
    VNewt = l_per_M**2 * r_per_M**-2 / 2 - r_per_M**-1
    plt.figure()
    plt.plot(r_per_M, Veff, 'r-', label=r'$V_{eff} \, \left( \frac{r}{M} \right)$')
    plt.plot(r_per_M, VNewt, 'b--', label=r'$V_\text{Newt} \, \left( \frac{r}{M} \right)$')
    plt.ylim([-0.04, 0.04])
    plt.xlabel(r'$\frac{r}{M}$')
    plt.ylabel(r'$V$')
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig('../latex/Figures/Veff.eps', format='eps')
    plt.show()

def plot_some_Veff(save=['yes','no']):

    lines = ['--', '-', ':', '-.']; kk = 0
    lines = ['-','-','-','-']
    l_per_M = np.array([8,  4.3, np.sqrt(12), 1])
    r_per_M = np.arange(2,50,0.1)
    plt.figure()

    for x in l_per_M:
        Veff = x**2 * r_per_M**-2 / 2 - r_per_M**-1 - x**2 * r_per_M**-3
        plt.plot(r_per_M, Veff, linestyle=lines[kk],  label=r'$V_{eff}, ~~ \frac{l}{M} = $'+str(x)[:4])
        kk += 1
    plt.ylim([-0.2, 0.4])
    plt.xlabel(r'$\frac{r}{M}$')
    plt.ylabel(r'$V$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    if save == 'yes': plt.savefig('../latex/Figures/Veff_tanti.eps', format='eps')
    plt.show()

''' Veff vs Newtonian V '''
# plot_Vs(4, 'yes')
plot_some_Veff('yes')
