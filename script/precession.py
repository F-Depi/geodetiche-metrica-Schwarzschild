# ::setlocal makeprg=cd\ script\ &&\ python\ precession.py
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.special import ellipk


SMALL_SIZE = 14
MEDIUM_SIZE = 15
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize


def plot_data_precession(foldername):

    filename = None
    for file in os.listdir(f'data/keep/{foldername}'):
        if file.endswith(".csv"):
            filename = file

    if not filename:
        print('No file found')
        exit()

    data = np.loadtxt(f'data/keep/{foldername}/{filename}', delimiter=',', skiprows=1)
    r = data[:, 1]
    phi = data[:, 2]

    plt.figure()
    plt.plot(phi, r, label='Precession')
    plt.title(f'Precession of massive particle in Schwarzschild metric\n{filename}')
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\frac{r}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.show()


def data_precession(l, E, h, other):

    filename = f'data/precession/l{l:.3f}_E{E:.5f}_h{h:.0e}{other}.csv'
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    r = data[:, 0]
    phi = data[:, 1]

    ## Clean the data from duplicates of the same turning point:
    ## sometimes the condition used in prec.c r_prev < r_turning < r_next is
    ## satisfied more than once in the same turn, so more than one point is
    ## saved as the turning point

    rmid = (np.max(r) + np.min(r)) / 2
    r_new = []
    phi_new = []
    for i in range(0, len(r) - 1):
        if r[i] > rmid and r[i+1] > rmid:
            continue
        if r[i] < rmid and r[i+1] < rmid:
            continue
        r_new.append(r[i])
        phi_new.append(phi[i])

    r = np.array(r_new)
    phi = np.array(phi_new)

    ## Precession with the outer turning points
    phi_outer = phi[r > rmid]
    prec_outer = []
    for i in range(1, len(phi_outer)):
        delta = (phi_outer[i] - phi_outer[i-1]) - 2 * np.pi
        prec_outer.append(delta)

    ## Precession with the inner turning points
    phi_inner = phi[r < rmid]
    prec_inner = []
    for i in range(1, len(phi_inner)):
        delta = (phi_inner[i] - phi_inner[i-1]) - 2 * np.pi
        prec_inner.append(delta)

    ## Precession from inner to outer and vice versa
    prec = []
    for i in range(1, len(phi)):
        delta = 2 * (phi[i] - phi[i-1] - np.pi)
        prec.append(delta)

    if r[0] < r[1]:
        prec_i2o = prec[::2]
        prec_o2i = prec[1::2]
    else:
        prec_i2o = prec[1::2]
        prec_o2i = prec[::2]

    return prec_outer, prec_inner, prec_i2o, prec_o2i


def analytic_precession(l, E):
    coefficients = [l**2, -l**2, 1, 2*E]
    roots = np.roots(coefficients)
    # We want the roots ordered in increasing order if they are all real
    # If there is just one real root, we want it to be the last element
    if not (np.all(np.isreal(roots))):
        print('There are complex roots')
        exit()

    roots = np.sort(roots)

    ## k^2 = (u2 - u1)/(u3 - u1)
    m = (roots[1] - roots[0])/(roots[2] - roots[0])
    K = ellipk(m)
    Delta_phi = 4 * K / np.sqrt(roots[2] - roots[0])
    return Delta_phi - 2 * np.pi


def plot_residuals(l, E, h, lloc=None, other=None, figname=None):

    # Get the data
    prec = analytic_precession(l, E)
    prec_outer, prec_inner, prec_i2o, prec_o2i = data_precession(l, E, h, other)

    # Residuals
    error_outer = prec_outer - prec 
    error_inner = prec_inner - prec 
    error_i2o = prec_i2o - prec
    error_o2i = prec_o2i - prec
    print(f'Error outer: {np.mean(error_outer):.3e} +/- {np.std(error_outer) / np.sqrt(len(error_outer)):.3e}')
    print(f'Error inner: {np.mean(error_inner):.3e} +/- {np.std(error_inner) / np.sqrt(len(error_inner)):.3e}')

    ## Normalize the residuals
    error_outer /= prec
    error_inner /= prec
    error_i2o /= prec
    error_o2i /= prec

    plt.figure()
    plt.plot(error_outer, linestyle='', marker='<', label=r'between $r_2$')
    plt.plot(error_inner, linestyle='', marker='.', label=r'between $r_1$')
    #plt.plot(error_i2o, linestyle='', marker='.', label=r'$r_1$ to $r_2$')
    #plt.plot(error_o2i, linestyle='', marker='.', label=r'$r_2$ to $r_1$')
    plt.axhline(0, color='black')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.xlabel('Revolutions')
    plt.ylabel('Relative error')
    plt.legend(loc=lloc)
    plt.title('Normalized precession residuals\n'rf'($\hat \ell = {l:.0f}$, $\mathcal{{E}} = {E:.3f}$, $h = {h:.0e}$)')
    plt.tight_layout()
    if figname is not None:
        plt.savefig(f'../latex/Figures/chapter2/{figname}')




## Write bisection to fine tune the energy
#for E in np.arange(-0.001, -0.01, -0.0005):
#    p = analytic_precession(2, E)
#    print(rf'E = {E:.4f}, $p = {p / np.pi:.3f} \pi$')



l = 5
E = -0.004
print(f'prec = {analytic_precession(l, E) / np.pi} pi')
plot_residuals(l, E, 1e-3, 'center right', other='_RK4')#, figname='prec1_res.eps')
plot_residuals(l, E, 1e-3, 'lower right', other='_RK4_corr')#, figname='prec1_res_corr.eps')


l = 3
E = -0.006
print(f'prec = {analytic_precession(l, E) / np.pi} pi')
plot_residuals(l, E, 1e-3,'center right', other='_RK4')#, figname='prec2_res.eps')
plot_residuals(l, E, 1e-3,'lower right', other='_RK4_corr')#, figname='prec2_res_corr.eps')


#plt.show()
