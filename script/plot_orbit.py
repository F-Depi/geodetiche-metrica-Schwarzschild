import matplotlib.pyplot as plt
import numpy as np
import sys
import os


''' DON'T MODIFY THIS FUNCTION '''
''' It's used by sim.sh to immediately plot the orbit '''

def plot_orbit1(l, E):
    
    if E == 'max':
        filename = f'l{l:.3f}_Emax.csv'
    elif E == 'min':
        filename = f'l{l:.3f}_Emin.csv'
    else:
        E = float(E)
        filename = f'l{l:.3f}_E{E:.5f}.csv'

    data = np.loadtxt(f'data/{filename}', delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]

    r_s = np.linspace(0, 2 * np.pi, 100)

    plt.figure()
    plt.plot(r * np.cos(phi), r * np.sin(phi),
             linestyle='-', marker='', markersize=0.5, label='orbit')
    plt.plot(r[0] * np.cos(phi[0]), r[0] * np.sin(phi[0]), 'ro', label='start')
    plt.plot(r[-1] * np.cos(phi[-1]), r[-1] * np.sin(phi[-1]), 'go', label='end')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label=r'$r_s$')
    plt.axis('equal')
    plt.title(rf'Massive Particle with $\hat \ell = {l}$, $\mathcal{{E}} = {E}$')
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.legend()
    plt.show()


def plot_orbit2(foldername, title):

    filename = None
    for file in os.listdir(f'data/keep/{foldername}'):
        if file.endswith(".csv"):
            filename = file

    if not filename:
        print('No file found')
        exit()

    data = np.loadtxt(f'data/keep/{foldername}/{filename}', delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]

    r_s = np.linspace(0, 2 * np.pi, 100)

    if title == '':
        l = filename[1:2]
        E = filename[8:15]
        title = rf'Massive Particle with $\hat \ell = {l}$, $\mathcal{{E}} = {E}$'

    plt.figure()
    plt.plot(r * np.cos(phi), r * np.sin(phi),
             linestyle='-', marker='.', markersize=1, label='orbit')
    plt.plot(r[0] * np.cos(phi[0]), r[0] * np.sin(phi[0]), 'ro', label='start')
    plt.plot(r[-1] * np.cos(phi[-1]), r[-1] * np.sin(phi[-1]), 'go', label='end')
    #plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label=r'$r_s$')
    plt.axis('equal')
    plt.title(title)
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.legend(loc='best')
    plt.show()


if len(sys.argv) == 3:
    l = float(sys.argv[1])
    E = sys.argv[2]
    plot_orbit1(l, E)


if len(sys.argv) == 2:
    foldername = sys.argv[1]
    plot_orbit2(foldername, '')

