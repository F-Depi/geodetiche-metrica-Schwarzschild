# ::setlocal makeprg=cd\ script\ &&\ python\ ch2_plots.py
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

def plot_potential():
    data = np.loadtxt('data/Veff.csv', delimiter=',', skiprows=1)
    r = data[:, 0]
    Veff = data[:, 1]

    plt.figure()
    plt.plot(r, Veff, label='Effective Potential')
    plt.show()


def plot_orbit():

    data = np.loadtxt('data/orbit.csv', delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    theta = data[:, 3]

    r_s = np.linspace(0, 2 * np.pi, 100)

    plt.figure()
    plt.plot(r * np.cos(phi), r * np.sin(phi),
             linestyle='', marker='.', label='orbit')
    plt.plot(r[0] * np.cos(phi[0]), r[0] * np.sin(phi[0]), 'ro', label='start')
    plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label='Event Horizon')
    plt.axis('equal')
    plt.title('Massive particle in Schwarzschild metric')
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.legend()
    plt.show()


#plot_potential()
plot_orbit()


