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


def plot_radial_infall(save=['yes','no']):

    # It's not just about plotting them, we also need to choose tau* and t* so
    # that the slopes appears together in the graph.
    # Also we have t(r) and not r(t)...
    # Get the start tau from r(tau)
    start_tau = 20
    tau = np.arange(0,start_tau,0.01)
    r_tau = (3/4)**(2/3) * (start_tau - tau)**(2/3)

    # Eval t(r), shift it so that it starts at t=0, plot (t(r), r)
    r = np.arange(2,r_tau[0],0.01)
    log_arg = (np.sqrt(r/2) + 1)/abs(np.sqrt(r/2) - 1)
    t = 2*( (-2/3) * (r/2)**(3/2) - 2*np.sqrt(r/2) + np.log(log_arg) )
    t = t - t[-1]

    # It wont' make any sense anyway because the functions where found from
    # r=infinity ad we are plotting a fall from a finite r.

    plt.figure()
    plt.plot(tau, r_tau, 'r-', label=r'$\frac{r}{M} \left( \frac{\tau}{M} \right)$')
    plt.plot(t, r, 'b-', label=r'$\frac{r}{M} \left( \frac{t}{M} \right)$')
    plt.xlabel(r'$\frac{time}{M}$')
    plt.ylabel(r'$\frac{r}{M}$')
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig('../latex/Figures/Veff.eps', format='eps')
    plt.show()

''' Veff vs Newtonian V '''
# plot_Vs(4, 'yes')
# plot_some_Veff('yes')

''' r(tau) vs r(t) '''
plot_radial_infall('no')
