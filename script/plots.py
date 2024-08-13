# ::setlocal makeprg=cd\ script\ &&\ python\ plots.py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar

SMALL_SIZE = 13
MEDIUM_SIZE = 18
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize


def fun_V_eff(r_per_M, l_per_M):
    return l_per_M**2 * r_per_M**-2 / 2 - r_per_M**-1 - l_per_M**2 * r_per_M**-3


def plot_Vs(l_per_M, save=['yes','no']):

    filename = '../latex/Figures/V_eff.eps'
    r_per_M = np.arange(2,40,0.1)
    Veff = fun_V_eff(r_per_M, l_per_M)
    label_eff = r'$V_{\rm eff} \, \left( \frac{r}{M} \right)$'
    VNewt = l_per_M**2 * r_per_M**-2 / 2 - r_per_M**-1
    label_Newt = r'$V_{\rm Newt} \, \left( \frac{r}{M} \right)$'

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


def plot_some_V_eff(save=['yes','no']):

    filename = '../latex/Figures/V_eff_tanti.eps'

    lines = ['--', '-', ':', '-.']; kk = 0
    lines = ['-','-','-','-']
    l_per_M = np.array([6,  4.3, np.sqrt(12), 1])
    l_str = ['6', '4.3', r'$\sqrt{12}$', '1']
    r_per_M = np.arange(2,50,0.1)
    plt.figure()

    for x in l_per_M:
        Veff = fun_V_eff(r_per_M, x)
        lab = r'$V_{\rm eff}, ~~ \frac{l}{M} = $' + l_str[kk]
        plt.plot(r_per_M, Veff, linestyle=lines[kk],  label=lab)
        kk += 1
    plt.xlim([0, 50])
    plt.ylim([-0.2, 0.36])
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
    plt.xlabel(r'$\frac{\rm time}{M}$')
    plt.ylabel(r'$\frac{r}{M}$', rotation=0)
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_V_eff_orbits(save=['yes','no']):

    filename = '../latex/Figures/V_eff_orbits.eps'
    right_lim = 70
    r_per_M = np.arange(2.1,right_lim,0.1)
    l_per_M = 4.1
    Veff = fun_V_eff(r_per_M, l_per_M)
    label_eff = r'$V_{\rm eff} \, (r)$'

    skrt = np.sqrt(1 - 12 / l_per_M**2)
    r_min = l_per_M**2 / 2 * (1 - skrt)
    V_1 = fun_V_eff(r_min, l_per_M)

    r_max = l_per_M**2 / 2 * (1 + skrt)
    V_2 = fun_V_eff(r_max, l_per_M)
    lab_e4 = r'$e_4 = V_{\rm eff} \, (r_{\rm max})$'

    e1 = 0.02
    x1_right =right_lim 
    lab_e1 = r'$e_1 =$'+str(e1)

    e2 = 0.005
    x2_left = root_scalar(lambda x: fun_V_eff(x, l_per_M) - e2, bracket=[r_min, r_max], method='bisect').root
    x2_right =right_lim 
    lab_e2 = r'$e_2 =$'+str(e2)

    e3 = -0.030
    x1_left = 0
    x3_left = root_scalar(lambda x: fun_V_eff(x, l_per_M) - e3, bracket=[r_min, r_max], method='bisect').root
    x3_right = root_scalar(lambda x: fun_V_eff(x, l_per_M) - e3, bracket=[r_max, right_lim], method='bisect').root
    lab_e3 = r'$e_3 =$'+str(e3)

    arrow_x = 40
    arrow_down = fun_V_eff(arrow_x, l_per_M)
    arrow_up = e2
    lab_arrow = r'$\left( \frac{\text{d} r}{\text{d} \tau} \right) ^2$  '

    plt.figure(figsize=(10,5))
    plt.plot(r_per_M, Veff, color='black', label=label_eff)
    plt.plot(r_min, V_1, 'rx', markersize=7)
    plt.hlines(e1, x1_left, x1_right, linestyle='--', color='b', label=lab_e1)
    plt.hlines(e2, x2_left, x2_right, linestyle='--', color='g', label=lab_e2)
    plt.hlines(e3, x3_left, x3_right, linestyle='--', color='orange', label=lab_e3)
    plt.plot(r_max, V_2, 'ro', label=lab_e4)
    plt.annotate('', xy=(arrow_x, arrow_down), xytext=(arrow_x, arrow_up),
             arrowprops=dict(arrowstyle='<->', lw=1.5, color='black'))
    plt.text(arrow_x, (arrow_up + arrow_down)/2, lab_arrow, fontsize=20, ha='right', va='center')
    plt.ylim([-0.06, 0.03])
    plt.xlim([0, right_lim])
    plt.xlabel(r'$\frac{r}{M}$')
    plt.ylabel(r'$V$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


''' Veff vs Newtonian V '''
#plot_Vs(4, 'yes')
#plot_some_V_eff('yes')

''' r(tau) vs r(t) '''
#plot_radial_infall('yes')

''' Veff vs values of e '''
plot_V_eff_orbits('yes')
