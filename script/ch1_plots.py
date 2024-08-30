# ::setlocal makeprg=cd\ script\ &&\ python\ ch1_plots.py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

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


def fun_W_eff(r_per_rs, M):
    return (1 - 1 / r_per_rs) / r_per_rs**2 / (2 * M)**2


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
        lab = r'$V_{\rm eff}, ~~ \frac{\ell}{M} = $' + l_str[kk]
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
    bottom_lim = -0.06

    r_per_M = np.arange(2.1,right_lim,0.01)
    l_per_M = 4.1
    Veff = fun_V_eff(r_per_M, l_per_M)
    label_eff = r'$V_{\rm eff} \, (r)$'

    skrt = np.sqrt(1 - 12 / l_per_M**2)
    r_min = l_per_M**2 / 2 * (1 - skrt)
    V_1 = fun_V_eff(r_min, l_per_M)

    # Different values of E
    r_max = l_per_M**2 / 2 * (1 + skrt)
    V_2 = fun_V_eff(r_max, l_per_M)
    lab_e4 = r'$\mathcal{E}_4 = V_{\rm eff} \, (r_+)$'

    e1 = 0.02
    x1_right =right_lim 
    lab_e1 = r'$\mathcal{E}_1 =$'+str(e1)

    e2 = 0.005
    x2_left = root_scalar(lambda x: fun_V_eff(x, l_per_M) - e2, bracket=[r_min, r_max], method='bisect').root
    x2_right =right_lim 
    lab_e2 = r'$\mathcal{E}_2 =$'+str(e2)

    e3 = -0.030
    x1_left = 0
    x3_left = root_scalar(lambda x: fun_V_eff(x, l_per_M) - e3, bracket=[r_min, r_max], method='bisect').root
    x3_right = root_scalar(lambda x: fun_V_eff(x, l_per_M) - e3, bracket=[r_max, right_lim], method='bisect').root
    lab_e3 = r'$\mathcal{E}_3 =$'+str(e3)

    # Vertical double arrow
    arrow_x = 40
    arrow_down = fun_V_eff(arrow_x, l_per_M)
    arrow_up = e2
    lab_arrow = r'$\left( \frac{\text{d} r}{\text{d} \tau} \right) ^2$  '


    ax = plt.figure(figsize=(10,5))
    plt.plot(r_per_M, Veff, color='black', label=label_eff)
    #plt.plot(r_min, V_1, 'rx', markersize=7)
    plt.vlines(r_min, bottom_lim, V_1, color='r', linestyle='--', linewidth=1)
    plt.vlines(r_max, bottom_lim, V_2, color='r', linestyle='--', linewidth=1)
    plt.hlines(e1, x1_left, x1_right, linestyle='--', color='b', label=lab_e1)
    plt.hlines(e2, x2_left, x2_right, linestyle='--', color='g', label=lab_e2)
    plt.hlines(e3, x3_left, x3_right, linestyle='--', color='orange', label=lab_e3)
    plt.plot(r_max, V_2, 'ro', label=lab_e4)
    plt.annotate('', xy=(arrow_x, arrow_down), xytext=(arrow_x, arrow_up),
             arrowprops=dict(arrowstyle='<->', lw=1.5, color='black'))
    plt.text(arrow_x, (arrow_up + arrow_down)/2, lab_arrow, fontsize=20, ha='right', va='center')

    plt.plot(x3_left, e3, color='orange', marker='o')
    plt.text(x3_left - 2, e3 - 5e-3, r'$P_1$', fontsize=15)

    plt.plot(x3_right, e3, color='orange', marker='o')
    plt.text(x3_right + 0.5, e3 - 5e-3, r'$P_2$', fontsize=15)

    plt.plot(x2_left, e2, color='green', marker='o')
    plt.text(x2_left + 0.5, e2 - 5e-3, r'$P_3$', fontsize=15)

    tiks = [0,10,20,30,40,50,60,70] + [r_min, r_max]
    lab_tiks = ['0','10','20','30','40','50','60','70'] + [r'$r_-$', r'$r_+$']
    plt.xticks(tiks, lab_tiks)

    plt.ylim([bottom_lim, 0.03])
    plt.xlim([0, right_lim])
    plt.xlabel(r'$\frac{r}{M}$')
    plt.ylabel(r'$V$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_W_eff(save=['yes','no']):

    filename = '../latex/Figures/W_eff.eps'
    right_lim = 15
    bottom_lim = -0.02
    upper_lim = 0.04
    M = 1
    plt.figure(figsize=(9, 4))
    r = np.arange(0.9, right_lim,0.001)
    Weff = M**2 * fun_W_eff(r, M)
    label_eff = r'$W_{\rm eff} \, (r)$'
    plt.plot(r, Weff, '-', linewidth=2, label=label_eff)
    plt.axhline(0, linewidth=2, color='black')
    plt.vlines(1, bottom_lim, 0, color='r', linestyle='--', linewidth=1)
    plt.vlines(1.5, bottom_lim, 1 / 27 / M**2, color='r', linestyle='--', linewidth=1)
    plt.xticks(list(plt.xticks()[0]) + [1])

    plt.ylim([bottom_lim, upper_lim])
    plt.xlim([0, right_lim])
    plt.xlabel(r'$\frac{r}{r_{\rm s}}$')
    plt.xlabel(r'$ r / r_s$')
    plt.ylabel(r'$M^2 ~ V$')
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_W_eff_tanti(save=['yes','no']):

    filename = '../latex/Figures/W_eff_tanti.eps'
    right_lim = 12
    bottom_lim = -0.02
    plt.figure()
    r = np.arange(0.9, right_lim, 0.001)
    lab_masses = ['4.43 mm', '1.48 km', '64,2 au']; kk = 0

    for M in [0.00443, 1480, 6.5e9 * 1480]:
        Weff = fun_W_eff(r, M)
        label_eff = r'$W_{\rm eff} \, (r),$ ' + lab_masses[kk]
        plt.plot(r, Weff, '-', linewidth=2, label=label_eff)
        kk += 1

    plt.axhline(0, color='black')
    plt.vlines(1, bottom_lim, 0, color='r', linestyle='--', linewidth=1)
    plt.xticks(list(plt.xticks()[0]) + [1])

    plt.xlim([0, right_lim])
    plt.yscale('log')
    plt.xlabel(r'$\frac{r}{r_{\rm s}}$')
    plt.ylabel(r'$V$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_W_eff_vs_b(save=['yes','no']):

    filename = '../latex/Figures/W_eff_vs_b.eps'
    right_lim = 15
    bottom_lim = 0
    upper_lim = 0.05
    M = 1

    b1 = 0.045
    x1 = 0
    lab_b1 = r'$M^2 / b^2 =$' + str(b1)

    b2 = 1 / 27
    x2 = 1.5
    lab_b2 = r'$M^2 / b^2 = 1 / 27$'

    b3 = 0.02
    x3 = root_scalar(lambda x: fun_W_eff(x, M) - b3,
                        bracket=[1.5, right_lim], method='bisect').root
    lab_b3 = r'$M^2 / b^2 =$' + str(b3)

    plt.figure(figsize=(10,5))
    r = np.arange(0.9, right_lim,0.001)
    Weff = M**2 * fun_W_eff(r, M)
    label_eff = r'$W_{\rm eff} \, (r)$'
    plt.plot(r, Weff, '-', linewidth=2, label=label_eff)

    # dashed lines and turning point
    plt.hlines(b1, x1, right_lim, linewidth=2, linestyle='--', color='orange', label=lab_b1)
    plt.hlines(b2, x2, right_lim, linewidth=2, linestyle='--', color='r', label=lab_b2)
    plt.hlines(b3, x3, right_lim, linewidth=2, linestyle='--', color='purple', label=lab_b3)
    plt.plot(x3, b3, color='purple', marker='o')
    plt.text(x3, b3 + 1e-3, r'$P$', fontsize=15)

    plt.ylim([bottom_lim, upper_lim])
    plt.xlim([0, right_lim])
    plt.xlabel(r'$\frac{r}{r_{\rm s}}$')
    plt.xlabel(r'$ r / r_s$')
    plt.ylabel(r'$M^2 ~ V$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def integrand_w(w, M_per_b):
    return (1 - w * w * (1 - 2 * w * M_per_b))**(-1/2)


def fun_phi_def(M_per_b):
    #to find w1 we want the integrand to go to 0. Rewriting it give a cubic
    # equation 2 M/b w^3 - w^2 + 1 = 0
    coefficients = [2 * M_per_b, -1, 0, 1]
    roots = np.roots(coefficients)
    real_roots = roots[np.isclose(roots.imag, 0)].real
    positive_roots = real_roots[real_roots > 0]
    # We want the bigger radius => the smaller w
    w1 = min(positive_roots)

    integral = quad(lambda w: integrand_w(w, M_per_b), 0, w1)[0]

    return 2 *integral - np.pi


def plot_light_deflection(save=['yes','no']):

    filename = '../latex/Figures/deflection_w.eps'

    limit = np.sqrt(1/27)
    lab_limit = r'$\frac{M}{b} = \frac{1}{\sqrt{27}}$'

    M_per_b = np.arange(0, limit, 0.0001)
    integral=[]
    for i in M_per_b:
        integral.append(fun_phi_def(i) / np.pi)
    lab_i = r'$\delta \phi$'

    plt.figure()
    plt.axvline(np.sqrt(1/27), color='r', linestyle='--', linewidth=2, label=lab_limit)
    plt.plot(M_per_b, integral, linewidth=2, label=lab_i)
    plt.xticks([0, 0.05, 0.1, 0.15, 0.2])

    plt.ylim([0, 2.5])
    plt.xlabel(r'$\frac{M}{b}$')
    plt.ylabel(r'$\frac{\delta \phi_{\rm def}}{\pi}$')
    plt.legend()
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


''' Veff vs Newtonian V '''
#plot_Vs(4, 'yes')
#plot_some_V_eff('yes')

''' r(tau) vs r(t) '''
#plot_radial_infall('yes')

''' Veff vs values of e '''
plot_V_eff_orbits('no')

''' Weff '''
#plot_W_eff('yes')
#plot_W_eff_tanti('no')
#plot_W_eff_vs_b('yes')

''' Light deflection '''
#plot_light_deflection('yes')


