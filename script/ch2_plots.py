# ::setlocal makeprg=cd\ script\ &&\ python\ ch2_plots.py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys

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


def plot_orbit(l, E):

    if l >= 0:
        filename = f'data/l{l:.3f}_E{E:.5f}.csv'
    else:
        filename = 'data/orbit.csv'

    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]

    r_s = np.linspace(0, 2 * np.pi, 100)

    plt.figure()
    plt.plot(r * np.cos(phi), r * np.sin(phi),
             linestyle='', marker='.', markersize=0.5, label='orbit')
    #plt.plot(r[0] * np.cos(phi[0]), r[0] * np.sin(phi[0]), 'ro', label='start')
    #plt.plot(r[-1] * np.cos(phi[-1]), r[-1] * np.sin(phi[-1]), 'go', label='end')
    plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label='Event Horizon')
    plt.axis('equal')
    plt.title('Massive particle in Schwarzschild metric')
    plt.xlabel(r'$\frac{x}{r_s}$')
    plt.ylabel(r'$\frac{y}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.legend()
    plt.show()


def precession(l, E):
    filename = f'data/l{l:.3f}_E{E:.5f}.csv'
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]

    plt.figure()
    plt.plot(phi, r, label='Precession')
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\frac{r}{r_s}$', rotation=0)
    plt.tight_layout()
    plt.show()


def animate_orbit(l, E):
    if l >= 0:
        filename = f'data/l{l:.3f}_E{E:.5f}.csv'
    else:
        filename = 'data/orbit.csv'

    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]

    r_s = np.linspace(0, 2 * np.pi, 100)

    # Calculate x and y coordinates
    x = r * np.cos(phi)
    y = r * np.sin(phi)

    # Create the figure and axis
    fig, ax = plt.subplots()

    # Initialize the plot with the first point of the trajectory
    scat = ax.scatter(x[0], y[0], c="b", s=10, label='Current Position')
    trajectory, = ax.plot([], [], c="r", lw=1, label='Trajectory')  # Line for trajectory
    plt.plot([0], [0], 'ko', label='black hole')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k--', label='Event Horizon')

    # Set axis labels
    ax.set(xlabel=r'$\frac{x}{r_s}$', ylabel=r'$\frac{y}{r_s}$')
    ax.legend()

    # Set the axis limits if needed
    ax.set_xlim(min(x) - 1, max(x) + 1)
    ax.set_ylim(min(y) - 1, max(y) + 1)

    # Update function for animation
    def update(frame):
        # Update the scatter plot with the current position
        scat.set_offsets([x[frame], y[frame]])
        
        # Update the trajectory line with all points up to the current frame
        trajectory.set_data(x[:frame+1], y[:frame+1])
        
        return scat, trajectory

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(tau), interval=1, blit=True)

    # Show the plot
    plt.show()

    filename = 'media/fall.mp4'
    ani.save(filename, dpi=400)



#plot_potential()
#l = float(sys.argv[1])
#E = float(sys.argv[2])
#plot_orbit(l, E)
#precession(l, E)
animate_orbit(-1, 0)





