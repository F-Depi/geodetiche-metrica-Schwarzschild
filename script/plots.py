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

    filename = '../latex/Figures/Veff.eps'
    r_per_M = np.arange(2,40,0.1)
    Veff = l_per_M**2 * r_per_M**-2 / 2 - r_per_M**-1 - l_per_M**2 * r_per_M**-3
    label_eff = r'$V_{eff} \, \left( \frac{r}{M} \right)$'
    VNewt = l_per_M**2 * r_per_M**-2 / 2 - r_per_M**-1
    label_Newt = r'$V_\text{Newt} \, \left( \frac{r}{M} \right)$'

    plt.figure()
    plt.plot(r_per_M, Veff, 'r-', label=label_eff)
    plt.plot(r_per_M, VNewt, 'b--', label=label_Newt)
    plt.ylim([-0.04, 0.04])
    plt.xlabel(r'$\frac{r}{M}$')
    plt.ylabel(r'$V$')
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_some_Veff(save=['yes','no']):

    filename = '../latex/Figures/Veff_tanti.eps'

    lines = ['--', '-', ':', '-.']; kk = 0
    lines = ['-','-','-','-']
    l_per_M = np.array([8,  4.3, np.sqrt(12), 1])
    r_per_M = np.arange(2,50,0.1)
    plt.figure()

    for x in l_per_M:
        Veff = x**2 * r_per_M**-2 / 2 - r_per_M**-1 - x**2 * r_per_M**-3
        lab = r'$V_{eff}, ~~ \frac{l}{M} = $'+str(x)[:4]
        plt.plot(r_per_M, Veff, linestyle=lines[kk],  label=lab)
        kk += 1
    plt.ylim([-0.2, 0.4])
    plt.xlabel(r'$\frac{r}{M}$')
    plt.ylabel(r'$V$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    if save == 'yes':
        plt.savefig(filename, format='eps')
    plt.show()


def plot_radial_infall(save=['yes','no']):

    filename = '../latex/Figures/radial_infall.eps'

    # It's not just about plotting them, we also need to choose tau* and t* so
    # that the slopes appears together in the graph.
    # Also we have t(r) and not r(t)...
    # Get the start tau from r(tau)
    start_tau =20
    tau = np.arange(0,start_tau + 0.01,0.01)
    r_tau = (9/2)**(1/3) * (start_tau - tau)**(2/3)
    label_tau = r'$(\tau, ~ r(\tau))$'

    # Eval t(r), shift it so that it starts at t=0, plot (t(r), r)
    r = np.arange(2.001,r_tau[0],0.0001)
    log_arg = (np.sqrt(r/2) + 1)/abs(np.sqrt(r/2) - 1)
    t = 2*( (-2/3) * (r/2)**(3/2) - 2*np.sqrt(r/2) + np.log(log_arg) )
    t += - t[-1]
    label_t = r'$(t(r), ~ r)$'

    # It wont' make any sense anyway because the functions where found from
    # r=infinity ad we are plotting a fall from a finite r.

    plt.figure()
    plt.axhline(2, color='black', linestyle='--')
    plt.plot(t, r, 'b-', label=label_t)
    plt.plot(tau, r_tau, 'r-', label=label_tau)
    plt.xlim([0,t[0]])
    plt.ylim([0, 12.5])
    plt.xlabel(r'$\frac{time}{M}$')
    plt.ylabel(r'$\frac{r}{M}$', rotation=0)
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()

''' Veff vs Newtonian V '''
# plot_Vs(4, 'yes')
# plot_some_Veff('yes')

''' r(tau) vs r(t) '''
plot_radial_infall('yes')
