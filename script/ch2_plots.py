# ::setlocal makeprg=cd\ script\ &&\ python\ ch2_plots.py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys
import os

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


def plot_orbit(foldername):

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

    plt.figure()
    plt.plot(r * np.cos(phi), r * np.sin(phi),
             linestyle='-', marker='.', markersize=1, label='orbit')
    plt.plot(r[0] * np.cos(phi[0]), r[0] * np.sin(phi[0]), 'ro', label='start')
    plt.plot(r[-1] * np.cos(phi[-1]), r[-1] * np.sin(phi[-1]), 'go', label='end')
    plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label='Event Horizon')
    plt.axis('equal')
    plt.title(f'Massive particle in Schwarzschild metric\n{filename}')
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.legend()
    plt.show()


def precession(foldername):

    filename = None
    for file in os.listdir(f'data/keep/{foldername}'):
        if file.endswith(".csv"):
            filename = f'data/keep/{foldername}/{file}'

    if not filename:
        print('No file found')
        exit()

    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]

    plt.figure()
    plt.plot(phi, r, label='Precession')
    plt.title(f'Precession of massive particle in Schwarzschild metric\n{filename}')
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\frac{r}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.show()


def animate_orbit(foldername):

    filename = None
    for file in os.listdir(f'data/keep/{foldername}'):
        if file.endswith(".csv"):
            filename = f'data/keep/{foldername}/{file}'

    if not filename:
        print('No file found')
        exit()

    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]

    r_s = np.linspace(0, 2 * np.pi, 100)

    # Calculate x and y coordinates
    x = r * np.cos(phi)
    y = r * np.sin(phi)

    min_ax = np.min([x, y]) * 1.1
    max_ax = np.max([x, y]) * 1.1
    print(min_ax, max_ax)
    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(9, 9))

    # Initialize the plot with the first point of the trajectory
    scat = ax.scatter(x[0], y[0], c="b", s=10, label='Particle')  # Point for current position
    trajectory, = ax.plot([], [], c="r", lw=1, label='Trajectory')  # Line for trajectory
    plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label='Event Horizon')

    pos_lab = rf'$\tau = {tau[0]:.2f}$, $t = {t[0]:.2f}$, $\Delta t = {t[0] - tau[0]:.2f}$'
    annotation = ax.text(0.05, 0.95, pos_lab, transform=ax.transAxes, fontsize=14,
                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set(xlabel=r'$\frac{x}{r_s}$', ylabel=r'$\frac{y}{r_s}$')
    #ax.set_xlim(np.min(x)-1, np.max(x)+1)
    plt.axis('equal')
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    ax.legend()

    # Set the axis limits if needed
    ax.set_xlim(min_ax, max_ax)
    ax.set_ylim(min_ax, max_ax)

    # Update function for animation
    def update(frame):
        # Update the scatter plot with the current position
        scat.set_offsets([x[frame], y[frame]])
        
        # Update the trajectory line with all points up to the current frame
        trajectory.set_data(x[:frame], y[:frame])

        pos_lab = rf'$\tau = {tau[frame]:.2f}$, $t = {t[frame]:.2f}$, $\Delta t = {t[frame] - tau[frame]:.2f}$'
        annotation.set_text(pos_lab)
        
        return scat, trajectory, annotation

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(tau), interval=2, blit=True, repeat=False)

    # Show the plot
    plt.show()

    filename = f'media/{foldername}.mp4'
    ani.save(filename, dpi=400, fps=60)



#plot_potential()


foldername = sys.argv[1]
plot_orbit(foldername)
#precession(foldername)
#animate_orbit(foldername)





