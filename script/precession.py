# ::setlocal makeprg=cd\ script\ &&\ python\ precession.py
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.special import ellipk


def data_precession(foldername):

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

    index_inners = []
    index_outers = []
    for i in range(1, len(r)-1):
        if r[i] > r[i-1] and r[i] > r[i+1]:
            index_inners.append(i)
        if r[i] < r[i-1] and r[i] < r[i+1]:
            index_outers.append(i)
    print(f'Data from: {filename}')

    for i in range(len(index_inners)-1):
        delta_phi = phi[index_inners[i+1]] - phi[index_inners[i]] - 2 * np.pi
        print(f'Data inners: {delta_phi:f}')

    for i in range(len(index_outers)-1):
        delta_phi = phi[index_outers[i+1]] - phi[index_outers[i]] - 2 * np.pi
        print(f'Data outers: {delta_phi:f}')

    plt.figure()
    plt.plot(phi, r, label='Precession')
    plt.title(f'Precession of massive particle in Schwarzschild metric\n{filename}')
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\frac{r}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.show()


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


l = 3
E = -0.006
prec = analytic_precession(l, E)
print(f'Analytical (l = {l:.3f}, E = {E:.5f}): {prec:f}')

data_precession('precession')
