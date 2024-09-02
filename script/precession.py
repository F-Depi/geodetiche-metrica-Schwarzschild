# ::setlocal makeprg=cd\ script\ &&\ python\ precession.py
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.special import ellipk


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


def data_precession(l, E, h):

    filename = f'data/precession/l{l:.3f}_E{E:.5f}_h{h:.0e}.csv'
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    r = data[:, 0]
    phi = data[:, 1]

    ## Precession with the outer turning points
    rmid = (r[0] + r[1]) / 2
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

    return prec_outer, prec_inner, prec


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


def plot_residuals(l, E, h):

    # Get the data
    prec = analytic_precession(l, E)
    prec_outer, prec_inner, prec_12 = data_precession(l, E, h)

    # Residuals
    error_outer = (prec_outer - prec) / prec
    error_inner = (prec_inner - prec) / prec
    error_12 = prec_12 - prec

    plt.figure()
    plt.plot(error_outer, linestyle='', marker='.', label='Outer')
    plt.plot(error_inner, linestyle='', marker='.', label='Inner')
    #plt.plot(error_12, label='Inner to outer')
    plt.axhline(0, color='black', label='Analytical')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.xlabel('Step')
    plt.ylabel('Precession')
    plt.legend()
    plt.title('Normalized precession residuals\n'rf'($\hat \ell = {l:.3f}$, $\mathcal{{E}} = {E:.5f}$, $h = {h:.0e}$)')
    plt.tight_layout()


## Write bisection to fine tune the energy
for E in np.arange(-0.001, -0.01, -0.0005):
    p = analytic_precession(2, E)
    print(rf'E = {E:.4f}, $p = {p / np.pi:.3f} \pi$')


#l = 3
#E = -0.006
#for h in [1e-1, 1e-2, 1e-3, 1e-4]:
#    plot_residuals(l, E, h)
#
#
#l = 3
#E = -0.01
#for h in [1e-1, 1e-2, 1e-3, 1e-4]:
#    plot_residuals(l, E, h)
#
#plt.show()

