import matplotlib.pyplot as plt
import numpy as np
import sys


''' DON'T MODIFY THIS FUNCTION '''
''' It's used by sim.sh to immediately plot the orbit '''

def plot_orbit(l, E):
    
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
    plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label='Event Horizon')
    plt.axis('equal')
    plt.title(rf'Massive Particle with $\hat \ell = {l}$, $\mathcal{{E}} = {E}$')
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.legend()
    plt.show()


l = float(sys.argv[1])
E = sys.argv[2]
plot_orbit(l, E)
