# ::setlocal makeprg=cd\ script\ &&\ python\ beautiful_plots.py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

SMALL_SIZE = 14
MEDIUM_SIZE = 15
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize

def plot_orbit(foldername, title, loc, figname):

    ## Find the file
    filename = None
    for file in os.listdir(f'data/keep/{foldername}'):
        if file.endswith(".csv"):
            filename = file
    if not filename:
        print('No file found')
        exit()

    ## Make the title
    if title == '':
        l = filename[1:4]
        E = filename[8:15]
        title = rf'Massive Particle with $\hat \ell \simeq {l}$, $\mathcal{{E}} \simeq {E}$'

    ## Prepare the legend location
    if loc == '': loc = 'best'

    data = np.loadtxt(f'data/keep/{foldername}/{filename}', delimiter=',', skiprows=1)
    r = data[:, 1]
    phi = data[:, 2]

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    r_s = np.linspace(0, 2 * np.pi, 100)

    plt.figure()
    plt.plot(x, y, linestyle='-', marker='', markersize=1, label='orbit')
    plt.plot(x[0], y[0], 'ro', label='start')
    plt.plot(x[-1], y[-1], 'go', label='end')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k-', label=r'$r = r_s$')
    plt.axis('equal')
    plt.title(title)
    plt.xlabel(r'$\hat x$')
    plt.ylabel(r'$\hat y$', rotation=0)
    plt.tight_layout()
    plt.legend(loc=loc)
    if figname != '': plt.savefig(figname)

def plot_zoom20(foldername, title, loc, figname):

    ## Find the file
    filename = None
    for file in os.listdir(f'data/keep/{foldername}'):
        if file.endswith(".csv"):
            filename = file
    if not filename:
        print('No file found')
        exit()

    ## Make the title
    if title == '':
        l = filename[1:4]
        E = filename[8:15]
        title = rf'Massive Particle with $\hat \ell \simeq {l}$, $\mathcal{{E}} \simeq {E}$'

    ## Prepare the legend location
    if loc == '': loc = 'best'

    data = np.loadtxt(f'data/keep/{foldername}/{filename}', delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    r_s = np.linspace(0, 2 * np.pi, 100)

    plt.figure()

    for i in np.arange(0, len(tau), 1):
        if abs(phi[i] % (2 * np.pi)) < 0.01 or abs(phi[i] % (2 * np.pi) - 2 * np.pi) < 0.01:
            plt.text(x[i], y[i], f'{tau[i]:.0f}', fontsize=8)

    plt.plot(x, y, linestyle='-', marker='', markersize=1, label='orbit')
    plt.plot(x[0], y[0], 'ro', label='start')
    plt.plot(x[-1], y[-1], 'go', label='end')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k-', label=r'$r = r_s$')
    #plt.xlim(5.335968 + 6e-6, 5.335968 + 9e-6)
    plt.title(title)
    plt.xlabel(r'$\hat x$')
    plt.ylabel(r'$\hat y$', rotation=0)
    plt.tight_layout()
    plt.legend(loc=loc)
    #if figname != '': plt.savefig(figname)


upfolder = 'fiori/'
savefolder = '../latex/Figures/appendixB/beautiful'
#plot_orbit(upfolder + '2', '', 'best', savefolder + '2.eps')
#plot_orbit(upfolder + '3', '', 'best', savefolder + '3.eps')
#plot_orbit(upfolder + '4', '', 'best', savefolder + '4.eps')
#plot_orbit(upfolder + '5', '', 'best', '')
#plot_orbit(upfolder + '10', '', 'best', savefolder + '10.eps')
#plot_orbit(upfolder + '20', '', 'best', savefolder + '20.eps')
#plot_orbit(upfolder + '20big', '', 'best', '')
plot_orbit('spirograph2', '', 'best', '')


plt.show()
