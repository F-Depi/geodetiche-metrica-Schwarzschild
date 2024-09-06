# ::setlocal makeprg=cd\ script\ &&\ python\ ch2_plots.py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import root_scalar
import numpy as np
import sys
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


def fun_V_eff(r, l):
    return (l**2 / r**2 - 1 / r - l**2 / r**3) / 2


def plot_scenario0(save=['yes','no']):

    filename = '../latex/Figures/chapter2/scenario0.eps'
    left_lim = 1
    right_lim = 300
    bottom_lim = -0.1
    up_lim = 0.015

    r = np.arange(1,right_lim,0.01)
    l = 1.2
    Veff = fun_V_eff(r, l)
    label_eff = r'$V_{\rm eff} \, (r)$'

    # Different values of E
    e1 = 0.005
    x1_left = left_lim
    x1_right = right_lim 
    lab_e1 = r'$\mathcal{E}_1 =$'+str(e1)

    e2 = -0.02
    x2_left = left_lim
    x2_right = root_scalar(lambda x: e2 - fun_V_eff(x, l), bracket=[1, right_lim], method='bisect').root
    lab_e2 = r'$\mathcal{E}_2 =$'+str(e2)

    e3 = -0.07
    x3_left = left_lim
    x3_right = root_scalar(lambda x: e3 - fun_V_eff(x, l), bracket=[1, right_lim], method='bisect').root
    lab_e3 = r'$\mathcal{E}_2 =$'+str(e3)

    ax = plt.figure(figsize=(10,4))
    plt.plot(r, Veff, color='black', label=label_eff)
    plt.hlines(e1, x1_left, x1_right, linestyle='--', color='b', label=lab_e1)
    plt.hlines(e2, x2_left, x2_right, linestyle='-', color='g', label=lab_e2)
    plt.hlines(e3, x3_left, x3_right, linestyle='--', color='g', label=lab_e3)

    plt.plot(x2_right, e2, color='green', marker='o')
    plt.text(x2_right - 0.4, e2 + 5e-3, r'$\hat r_2  $', fontsize=15)

    plt.plot(x3_right, e3, color='green', marker='o')
    plt.text(x3_right - 0.4, e3 + 5e-3, r'$\hat r_2  $', fontsize=15)

    plt.xscale('log')
    tiks = [1,10,30,100,300]
    lab_tiks = ['1','10','30','100','300']
    plt.xticks(tiks, lab_tiks)

    plt.ylim([bottom_lim, up_lim])
    plt.xlim([left_lim, right_lim])
    plt.title(rf'Different orbit for $\hat \ell = {l:.1f}$')
    plt.xlabel(r'$\hat r$')
    plt.ylabel(r'$V$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_scenario1(save=['yes','no']):

    filename = '../latex/Figures/chapter2/scenario1.eps'
    left_lim = 1
    right_lim = 300
    bottom_lim = -0.06
    up_lim = 0.015

    r = np.arange(1,right_lim,0.01)
    l = 1.9
    Veff = fun_V_eff(r, l)
    label_eff = r'$V_{\rm eff} \, (r)$'

    skrt = np.sqrt(1 - 3 / l**2)
    r_max = l**2 * (1 - skrt)
    V_max = fun_V_eff(r_max, l)

    # Different values of E
    r_min = l**2 * (1 + skrt)
    V_min = fun_V_eff(r_min, l)

    e1 = 0.005
    x1_left = left_lim
    x1_right = right_lim 
    lab_e1 = r'$\mathcal{E}_1 =$'+str(e1)

    e2 = -0.015
    x2_left = left_lim
    x2_right = root_scalar(lambda x: e2 - fun_V_eff(x, l), bracket=[r_min, right_lim], method='bisect').root
    lab_e2 = r'$\mathcal{E}_2 =$'+str(e2)

    e3 = -0.030
    x3_left = root_scalar(lambda x: fun_V_eff(x, l) - e3, bracket=[r_max, r_min], method='bisect').root
    x3_right = root_scalar(lambda x: fun_V_eff(x, l) - e3, bracket=[r_min, right_lim], method='bisect').root
    lab_e3 = r'$\mathcal{E}_3 =$'+str(e3)

    e4 = -0.035
    x4_left = left_lim
    x4_right = root_scalar(lambda x: e4 - fun_V_eff(x, l), bracket=[1, r_max], method='bisect').root
    lab_e4 = r'$\mathcal{E}_4 =$'+str(e4)

    e5 = -0.053
    x5_left = left_lim
    x5_right = root_scalar(lambda x: e5 - fun_V_eff(x, l), bracket=[1, r_max], method='bisect').root
    lab_e5 = r'$\mathcal{E}_5 =$'+str(e5)

    ax = plt.figure(figsize=(10,5))
    plt.plot(r, Veff, color='black', label=label_eff)
    plt.vlines(r_max, bottom_lim, V_max, color='r', linestyle='--', linewidth=1)
    plt.vlines(r_min, bottom_lim, V_min, color='r', linestyle='--', linewidth=1)
    plt.hlines(e1, x1_left, x1_right, linestyle='--', color='b', label=lab_e1)
    plt.hlines(e2, x2_left, x2_right, linestyle='-', color='g', label=lab_e2)
    plt.hlines(e3, x3_left, x3_right, linestyle='--', color='purple', label=lab_e3)
    plt.hlines(e4, x4_left, x4_right, linestyle='--', color='g', label=lab_e4)
    plt.hlines(e5, x5_left, x5_right, linestyle=':', color='g', label=lab_e5)

    plt.plot(x2_right, e2, color='green', marker='o')
    plt.text(x2_right, e2 + 2e-3, r'$\hat r_2  $', fontsize=15)

    plt.plot(x3_left, e3, color='purple', marker='o')
    plt.text(x3_left, e3 + 2e-3, r'$\hat r_1$', fontsize=15)

    plt.plot(x3_right, e3, color='purple', marker='o')
    plt.text(x3_right - 0.6, e3 + 2e-3, r'$\hat r_2$', fontsize=15)

    plt.plot(x4_right, e4, color='g', marker='o')
    plt.text(x4_right - 0.2, e4 + 2e-3, r'$\hat r_2  $', fontsize=15)

    plt.plot(x5_right, e5, color='g', marker='o')
    plt.text(x5_right - 0.2, e5 + 2e-3, r'$\hat r_2  $', fontsize=15)

    plt.xscale('log')
    tiks = [1,10,30,100,300] + [r_max, r_min]
    lab_tiks = ['1','10','30','100','300'] + [r'$\hat r_{\rm max}$', r'$\hat r_{\rm min}$']
    plt.xticks(tiks, lab_tiks)

    plt.ylim([bottom_lim, up_lim])
    plt.xlim([left_lim, right_lim])
    plt.title(rf'Different orbit for $\hat \ell = {l:.1f}$')
    plt.xlabel(r'$\hat r$')
    plt.ylabel(r'$V$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_scenario2(save=['yes','no']):

    filename = '../latex/Figures/chapter2/scenario2.eps'
    left_lim = 1
    right_lim = 300
    bottom_lim = -0.06
    up_lim = 0.08

    r = np.arange(1,right_lim,0.01)
    l = 2.15
    Veff = fun_V_eff(r, l)
    label_eff = r'$V_{\rm eff} \, (r)$'

    skrt = np.sqrt(1 - 3 / l**2)
    r_max = l**2 * (1 - skrt)
    V_max = fun_V_eff(r_max, l)

    # Different values of E
    r_min = l**2 * (1 + skrt)
    V_min = fun_V_eff(r_min, l)
    lab_e4 = r'$\mathcal{E}_4 = V_{\rm eff} \, (\hat r_{\rm min})$'

    e1 = 0.065
    x1_left = left_lim
    x1_right =right_lim 
    lab_e1 = r'$\mathcal{E}_1 =$'+str(e1)

    e2 = 0.025
    x2_left = root_scalar(lambda x: e2 - fun_V_eff(x, l), bracket=[r_max, r_min], method='bisect').root
    x2_right =right_lim 
    lab_e2 = r'$\mathcal{E}_2 =$'+str(e2)

    e3 = -0.02
    x3_left = root_scalar(lambda x: e3 - fun_V_eff(x, l), bracket=[r_max, r_min], method='bisect').root
    x3_right = root_scalar(lambda x: e3 - fun_V_eff(x, l), bracket=[r_min, right_lim], method='bisect').root
    lab_e3 = r'$\mathcal{E}_3 =$'+str(e3)

    e4 = +0.01
    x4_left = left_lim
    x4_right = root_scalar(lambda x: e4 - fun_V_eff(x, l), bracket=[1, r_max], method='bisect').root
    lab_e4 = r'$\mathcal{E}_4 =$'+str(e4)

    e5 = -0.04
    x5_left = left_lim
    x5_right = root_scalar(lambda x: e5 - fun_V_eff(x, l), bracket=[1, r_max], method='bisect').root
    lab_e5 = r'$\mathcal{E}_4 =$'+str(e5)

    ax = plt.figure(figsize=(10,5))
    plt.plot(r, Veff, color='black', label=label_eff)
    plt.vlines(r_max, bottom_lim, V_max, color='r', linestyle='--', linewidth=1)
    plt.vlines(r_min, bottom_lim, V_min, color='r', linestyle='--', linewidth=1)
    plt.hlines(e1, x1_left, x1_right, linestyle='--', color='b', label=lab_e1)
    plt.hlines(e2, x2_left, x2_right, linestyle='--', color='orange', label=lab_e2)
    plt.hlines(e3, x3_left, x3_right, linestyle='--', color='purple', label=lab_e3)
    plt.hlines(e4, x4_left, x4_right, linestyle='-', color='green', label=lab_e4)
    plt.hlines(e5, x5_left, x5_right, linestyle='--', color='green', label=lab_e5)

    plt.plot(x2_left, e2, color='orange', marker='o')
    plt.text(x2_left, e2 + 3e-3, r'$\hat r_1  $', fontsize=15)

    plt.plot(x3_left, e3, color='purple', marker='o')
    plt.text(x3_left , e3 + 3e-3, r'$\hat r_1$', fontsize=15)

    plt.plot(x3_right, e3, color='purple', marker='o')
    plt.text(x3_right, e3 + 3e-3, r'$\hat r_2$', fontsize=15)

    plt.plot(x4_right, e4, color='green', marker='o')
    plt.text(x4_right + 0.05, e4 + 3e-3, r'$\hat r_2  $', fontsize=15)

    plt.plot(x5_right, e5, color='green', marker='o')
    plt.text(x5_right + 0.05, e5 + 3e-3, r'$\hat r_2  $', fontsize=15)

    plt.xscale('log')
    tiks = [1,10,30,100,300] + [r_max, r_min]
    lab_tiks = ['1','10','30','100','300'] + [r'$\hat r_{\rm max}$', r'$\hat r_{\rm min}$']
    plt.xticks(tiks, lab_tiks)

    plt.ylim([bottom_lim, up_lim])
    plt.xlim([left_lim, right_lim])
    plt.title(rf'Different orbits for $\hat \ell = {l:.1f}$')
    plt.xlabel(r'$\hat r$')
    plt.ylabel(r'$V$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save == 'yes': plt.savefig(filename, format='eps')
    plt.show()


def plot_potential():
    data = np.loadtxt('data/Veff.csv', delimiter=',', skiprows=1)
    r = data[:, 0]
    Veff = data[:, 1]

    plt.figure()
    plt.plot(r, Veff, label='Effective Potential')
    plt.show()


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
        title = rf'Massive Particle with $\hat \ell = {l}$, $\mathcal{{E}} = {E}$'

    ## Prepare the legend location
    if loc == '': loc = 'best'

    data = np.loadtxt(f'data/keep/{foldername}/{filename}', delimiter=',', skiprows=1)
    tau = data[:, 0]
    r = data[:, 1]
    phi = data[:, 2]
    t = data[:, 3]
    tau = tau[r > 1]
    phi = phi[r > 1]
    t = t[r > 1]
    r = r[r > 1]

    r_s = np.linspace(0, 2 * np.pi, 100)

    plt.figure()
    plt.plot(r * np.cos(phi), r * np.sin(phi),
             linestyle='-', marker='', markersize=1, label='orbit')
    plt.plot(r[0] * np.cos(phi[0]), r[0] * np.sin(phi[0]), 'ro', label='start')
    plt.plot(r[-1] * np.cos(phi[-1]), r[-1] * np.sin(phi[-1]), 'go', label='end')
    plt.plot(np.cos(r_s), np.sin(r_s), 'k-', label=r'$r = r_s$')
    plt.axis('equal')
    plt.title(title)
    plt.xlabel(r'$\hat x$')
    plt.ylabel(r'$\hat y$', rotation=0)
    plt.tight_layout()
    plt.legend(loc=loc)
    if figname != '': plt.savefig(figname)


def animate_orbit(foldername):

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

    # Calculate x and y coordinates
    x = r * np.cos(phi)
    y = r * np.sin(phi)

    min_ax = np.min([x, y]) * 1.1
    max_ax = np.max([x, y]) * 1.1

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

    filename = f'media/{foldername}2.mp4'
    #ani.save(filename, writer='ffmpeg', fps=90, dpi=200, extra_args=['-vcodec', 'libx264'])


def fun_tau_r(r):
    return - 2 / 3 * r**(3/2)


def fun_t_r(r):
    return - 2 / 3 * r**(3/2) - 2 * r**(1/2) + np.log(np.abs( (r**(1/2) + 1) / (r**(1/2) - 1) ))


def plt_tvstau(foldername, h, save=['yes','no']):

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
    t = data[:, 3]

    tau = tau[r > 1]
    t = t[r > 1]
    r = r[r > 1]

    tau_int_const = tau[0] - fun_tau_r(r[0])
    analytic_tau = tau_int_const + fun_tau_r(r)

    t_int_const = t[0] - fun_t_r(r[0])
    analytic_t = t_int_const + fun_t_r(r)


    plt.figure()
    plt.plot(r, t, marker='.', linestyle='', label=r'$\hat t$ vs $\hat r$')
    plt.plot(r, analytic_t, label=r'$\hat t_{\rm analytic}$($\hat r)$')
    plt.plot(r, tau, marker='.', linestyle='', label=r'$\hat \tau$ vs $\hat r$')
    plt.plot(r, analytic_tau, label=r'$\hat \tau_{\rm analytic}$($\hat r)$')
    plt.title(f'Radial Infall')
    plt.xlabel(r'$\hat r$')
    plt.ylabel('time')
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.legend()
    if save == 'yes':
        plt.savefig(f'../latex/Figures/chapter2/radial_infall.eps')

    r = r[tau > 0.1]
    t = t[tau > 0.1]
    analytic_t = analytic_t[tau > 0.1]
    analytic_tau = analytic_tau[tau > 0.1]
    tau = tau[tau > 0.1]

    #t_res = abs(t / analytic_t - 1)
    #tau_res = abs(tau / analytic_tau - 1)
    t_res = abs(t - analytic_t)
    tau_res = abs(tau - analytic_tau)


    plt.figure(figsize=(7,5))
    plt.plot(r, t_res, linestyle='', marker='.',
             label=r'$\hat t - \hat t_{\rm analytic}$')
    plt.yscale('log')
    plt.xlabel(r'$\hat r$')
    plt.ylabel('Residuals')
    plt.title(rf'Normalized Residual of $\hat t$, $h = {h:.1e}$')
    plt.tight_layout()
    plt.legend()
    if save == 'yes':
        plt.savefig('../latex/Figures/chapter2/t_res.png')

    plt.figure(figsize=(7,5))
    plt.plot(r, tau_res, linestyle='', marker='.',
             label=r'$\hat \tau - \hat \tau_{\rm analytic}$')
    plt.xlabel(r'$\hat r$')
    plt.ylabel('Residuals')
    plt.title(rf'Normalized Residual of $\hat \tau$, $h = {h:.1e}$')
    plt.tight_layout()
    plt.legend()
    if save == 'yes':
        plt.savefig('../latex/Figures/chapter2/tau_res.png')


def plt_tvstau_multi(hs, save=['yes','no']):

    fig1, ax1 = plt.subplots(figsize=(7,5))
    fig2, ax2 = plt.subplots(figsize=(7,5))

    for h in hs:

        filename = f'data/keep/rad_infall_stability/l0_E0_h{h:.0e}.csv'
        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        tau = data[:, 0]
        r = data[:, 1]
        t = data[:, 3]

        tau = tau[r > 1]
        t = t[r > 1]
        r = r[r > 1]

        tau_int_const = tau[0] - fun_tau_r(r[0])
        analytic_tau = tau_int_const + fun_tau_r(r)

        t_int_const = t[0] - fun_t_r(r[0])
        analytic_t = t_int_const + fun_t_r(r)

        r = r[tau > 0.1]
        t = t[tau > 0.1]
        analytic_t = analytic_t[tau > 0.1]
        analytic_tau = analytic_tau[tau > 0.1]
        tau = tau[tau > 0.1]

        #t_res = abs(t / analytic_t - 1)
        #tau_res = abs(tau / analytic_tau - 1)
        t_res = abs(t - analytic_t) / t
        tau_res = (tau - analytic_tau) / tau


        h = float(h)
        ax1.plot(r, t_res, linestyle='', marker='.',
                 label=rf'$h = {h:.0e}$')

        ax2.plot(r, tau_res, linestyle='', marker='.',
                 label=rf'$h = {h:.0e}$')

    ax1.set_yscale('log')
    ax1.set_xlabel(r'$\hat r$')
    ax1.set_ylabel('Residuals')
    ax1.set_title(r'Normalized Residual: $\hat t / \hat t_{\rm analytic} - 1$')
    fig1.tight_layout()
    ax1.legend(loc='upper right')
    if save == 'yes':
        fig1.savefig('../latex/Figures/chapter2/t_res_multi.png')

    ax2.set_xlabel(r'$\hat r$')
    ax2.set_ylabel('Residuals')
    ax2.set_title(r'Normalized Residual: $\hat \tau / \hat \tau_{\rm analytic} - 1$')
    fig2.tight_layout()
    ax2.legend(loc='upper right')
    if save == 'yes':
        fig2.savefig('../latex/Figures/chapter2/tau_res_multi.png')


def check_stability_infall2(hs, save=['yes','no']):
    foldername = 'data/keep/infall2_stability/l1.8_E-0.042_h'
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create an inset axes for the zoomed-in region
    ax_inset1 = inset_axes(ax, width="30%", height="30%", loc='upper right')
    ax_inset2 = inset_axes(ax, width="30%", height="30%", loc='lower right')

    colors = ["#00202e","#003f5c","#2c4875","#8a508f","#bc5090","#ff6361","#ff8531","#ffa600","#ffd380"]
    colors = ["#ffd380","#ffa600","#ff8531","#ff6361","#bc5090","#8a508f","#2c4875","#003f5c","#00202e"]

    tau = []
    x = []
    y = []

    for h, col in zip(hs, colors):
        # Construct the filename
        filename = foldername + f'{h:.1e}.csv'

        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        tau = data[:, 0]
        r = data[:, 1]
        phi = data[:, 2]

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        # Plot the main data
        ax.plot(x, y, color=col, label=f'$h = {h:.1e}$')
        # Plot the zoomed-in data
        ax_inset1.plot(x, y, color=col)
        ax_inset2.plot(x, y, color=col)


    ## ZOOM 1
    # Define the region to zoom in
    x1 = - 4.044 - 0.00095
    x2 = - 4.044 - 0.00075 
    y1 =   5.477 - 0.00010
    y2 =   5.477 + 0.00000 
    ax_inset1.set_xlim(x1, x2)
    ax_inset1.set_ylim(y1, y2)

    ax_inset1.set_xticks([-9*1e-4 - 4.044, -8*1e-4 - 4.044])
    ax_inset1.set_yticklabels([])
    ax_inset1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Add a rectangle to highlight the zoomed-in area in the main plot
    #rect = plt.Rectangle((x1, y1), x2-x1, y2-y1, edgecolor='red', facecolor='none', linestyle='--')
    #ax.add_patch(rect)

    # Mark the inset on the main plot
    mark_inset(ax, ax_inset1, loc1=2, loc2=3, fc="none", ec="0.5")

    ## ZOOM 2
    #x1 = -1.976 - 51.25e-5
    #x2 = -1.976 - 49.75e-5
    #y1 = 0.2359 + 24e-6
    #y2 = 0.2359 + 32e-6
    x1 = 0.1834 + 0e-5
    x2 = 0.1834 + 5e-5
    y1 = 2.3002 + 5.7e-5
    y2 = 2.3002 + 8.3e-5
    ax_inset2.set_xlim(x1, x2)
    ax_inset2.set_ylim(y1, y2)

    N = 1932
    x_ann = 0.1834323
    y_ann = 2.3002672
    ax_inset2.annotate(rf'$\hat \tau \simeq {tau[N]:.2f}$', xy=(x_ann, y_ann),
                       xytext=(x_ann - 0.5e-5, y_ann - 1e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    N = 70146
    x_ann = 0.18341669
    y_ann = 2.3002783
    lab = rf'$\hat \tau \simeq {tau[N]:.2f}$'
    ax_inset2.annotate(lab, xy=(x_ann + 0.8e-5, y_ann),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    ax_inset2.annotate('               ', xy=(x_ann + 0.8e-5, y_ann - 0.7e-5),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    ax_inset2.annotate('                 ', xy=(x_ann + 0.8e-5, y_ann - 0.38e-5),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    ax_inset2.annotate('               ', xy=(x_ann + 0.8e-5, y_ann - 1.7e-5),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)

    ax_inset2.set_xticks([1e-5 + 0.1834, 2e-5 + 0.1834, 3e-5 + 0.1834, 4e-5 + 0.1834])
    ax_inset2.set_yticklabels([])
    ax_inset2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax_inset2.xaxis.set_ticks_position('top') 
    #ax_inset2.axis('equal')

    mark_inset(ax, ax_inset2, loc1=2, loc2=3, fc="none", ec="0.5")

    ## Main plot
    r_s = np.linspace(0, 2 * np.pi, 100)
    ax.plot(np.cos(r_s), np.sin(r_s), 'k-', label=r'$\hat r = \hat r_s$')
    ax.plot(x[0], y[0], 'ro', label='Start')
    ax.plot(x[-1], y[-1], 'go', label='End')

    ax.axis('equal')
    ax.set_xlim(-12.5, 15)
    ax.set_title(r'Stability of the numerical solution, $\hat r$ and $\phi$')
    ax.set_xlabel(r'$\hat x$')
    ax.set_ylabel(r'$\hat y$', rotation=0)
    ax.legend(loc='upper left')
    plt.tight_layout()

    if save == 'yes':
        plt.savefig(f'../latex/Figures/chapter2/stability_infall2.eps')

    plt.show()


def check_stability_volevi(hs, save=['yes','no']):
    foldername = 'data/keep/volevi_stability/l1.9_E-0.02367591_h'
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create an inset axes for the zoomed-in region
    ax_inset1 = inset_axes(ax, width="30%", height="30%", loc='upper right')
    ax_inset2 = inset_axes(ax, width="30%", height="30%", loc='lower right')

    colors = ["#00202e","#003f5c","#2c4875","#8a508f","#bc5090","#ff6361","#ff8531","#ffa600","#ffd380"]
    colors = ["#ffd380","#ffa600","#ff8531","#ff6361","#bc5090","#8a508f","#2c4875","#003f5c","#00202e"]

    tau = []
    x = []
    y = []

    for h, col in zip(hs, colors):
        # Construct the filename
        filename = foldername + f'{h:.1e}.csv'

        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        tau = data[:, 0]
        r = data[:, 1]
        phi = data[:, 2]

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        # Plot the main data
        ax.plot(x, y, color=col, label=f'$h = {h:.1e}$')
        # Plot the zoomed-in data
        ax_inset1.plot(x, y, color=col)
        ax_inset2.plot(x, y, color=col)


    ## ZOOM 1
    # Define the region to zoom in
    x1 = - 4.044 - 0.00095
    x2 = - 4.044 - 0.00075 
    y1 =   5.477 - 0.00010
    y2 =   5.477 + 0.00000 
    ax_inset1.set_xlim(x1, x2)
    ax_inset1.set_ylim(y1, y2)

    ax_inset1.set_xticks([-9*1e-4 - 4.044, -8*1e-4 - 4.044])
    ax_inset1.set_yticklabels([])
    ax_inset1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Add a rectangle to highlight the zoomed-in area in the main plot
    #rect = plt.Rectangle((x1, y1), x2-x1, y2-y1, edgecolor='red', facecolor='none', linestyle='--')
    #ax.add_patch(rect)

    # Mark the inset on the main plot
    mark_inset(ax, ax_inset1, loc1=2, loc2=3, fc="none", ec="0.5")

    ## ZOOM 2
    #x1 = -1.976 - 51.25e-5
    #x2 = -1.976 - 49.75e-5
    #y1 = 0.2359 + 24e-6
    #y2 = 0.2359 + 32e-6
    x1 = 0.1834 + 0e-5
    x2 = 0.1834 + 5e-5
    y1 = 2.3002 + 5.7e-5
    y2 = 2.3002 + 8.3e-5
    ax_inset2.set_xlim(x1, x2)
    ax_inset2.set_ylim(y1, y2)

    N = 1932
    x_ann = 0.1834323
    y_ann = 2.3002672
    ax_inset2.annotate(rf'$\hat \tau \simeq {tau[N]:.2f}$', xy=(x_ann, y_ann),
                       xytext=(x_ann - 0.5e-5, y_ann - 1e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    N = 70146
    x_ann = 0.18341669
    y_ann = 2.3002783
    lab = rf'$\hat \tau \simeq {tau[N]:.2f}$'
    ax_inset2.annotate(lab, xy=(x_ann + 0.8e-5, y_ann),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    ax_inset2.annotate('               ', xy=(x_ann + 0.8e-5, y_ann - 0.7e-5),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)
    ax_inset2.annotate('               ', xy=(x_ann + 0.8e-5, y_ann - 1.7e-5),
                       xytext=(x_ann - 1.5e-5, y_ann - 1.5e-5),
                       arrowprops=dict(facecolor='black', arrowstyle='-'),)

    ax_inset2.set_xticks([1e-5 + 0.1834, 2e-5 + 0.1834, 3e-5 + 0.1834, 4e-5 + 0.1834])
    ax_inset2.set_yticklabels([])
    ax_inset2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax_inset2.xaxis.set_ticks_position('top') 
    #ax_inset2.axis('equal')

    mark_inset(ax, ax_inset2, loc1=2, loc2=3, fc="none", ec="0.5")

    ## Main plot
    r_s = np.linspace(0, 2 * np.pi, 100)
    ax.plot(np.cos(r_s), np.sin(r_s), 'k-', label=r'$\hat r = \hat r_s$')
    ax.plot(x[0], y[0], 'ro', label='Start')
    ax.plot(x[-1], y[-1], 'go', label='End')

    ax.axis('equal')
    ax.set_xlim(-12, 29)
    ax.set_title(r'Stability of the numerical solution, $\hat r$ and $\phi$')
    ax.set_xlabel(r'$\hat x$')
    ax.set_ylabel(r'$\hat y$', rotation=0)
    ax.legend(loc='upper left')
    plt.tight_layout()

    if save == 'yes':
        plt.savefig(f'../latex/Figures/chapter2/stability_volevi.eps')

    plt.show()


def check_time_infall2(hs, save=['yes','no']):

    foldername = 'data/keep/infall2_stability/l1.8_E-0.042_h'
    colors = ["#ffd380","#ffa600","#ff8531","#ff6361","#bc5090","#8a508f","#2c4875","#003f5c","#00202e"]

    fig, ax = plt.subplots()

    for h, col in zip(hs, colors):

        filename = foldername + f'{h:.1e}.csv'
        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        tau = data[:, 0]
        r = data[:, 1]
        t = data[:, 3]

        # Plot the main data
        ax.plot(tau, r, color=col, label=rf'$\hat t (\hat r), ~ h = {h:.1e}$')
        #ax.plot(tau, t, color=col, label=rf'$\hat t (\hat r), ~ h = {h:.1e}$')
        print(f'{h:.1e}: {t[-1] - tau[-1]:.3f}')
    

    ax.set_title(r'Stability of the numerical solution, $\hat t$ and $\hat \tau$')
    ax.set_xlabel(r'$\hat \tau$')
    ax.set_ylabel(r'$\hat r, ~ \hat \tau$', rotation=0)
    #ax.legend()
    plt.tight_layout()

    # Save the figure if needed
    save = 'yes'
    if save == 'yes':
        plt.savefig(f'../latex/Figures/chapter2/t_stability_infall2.eps')

    plt.show()


def check_circular(foldername, h, figname):
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

    N = 1
    omega = [(phi[i+N] - phi[i]) / (t[i+N] - t[i]) for i in range(len(tau)-N)]
    omega = np.array(omega)

    ## Omega - Omega_analytic / Omega_analytic = sqrt(2 r^3) Omega - 1
    omega_analytic = (1 / 2 / r[0]**3)**(1/2)
    print(omega_analytic)
    const = np.sqrt(2 * r[0]) * r[0] 
    err = const * omega - 1

    lab = r'$\frac{\Omega - \Omega_{\rm analytic}}{\Omega_{\rm analytic}}$'
    plt.figure()
    plt.plot(t[N:], err, 'b ', marker='.', label=lab)
    plt.title(rf'Normalized residual for $\Omega$, $\hat \ell = 5$')
    plt.xlabel(r'$\hat t$')
    plt.tight_layout()
    if figname != '': plt.savefig(figname, dpi=200)
    plt.legend(fontsize=18, loc='lower left')



''' Potential check '''
#plot_potential()


''' Different orbits '''
#plot_scenario0('yes')
#plot_scenario1('yes')
#plot_scenario2('yes')


fold = '../latex/Figures/chapter2/'

''' radial infall '''
#plt_tvstau('rad_infall2e-3_RK4', 2e-3, 'no')
#plt_tvstau('rad_infall1e-3_RK4', 1e-3, 'yes')
#plt_tvstau('rad_infall5e-4_RK4', 5e-4, 'no')
plt_tvstau_multi([2e-3, 1e-3, 4e-4, 2e-4, 1e-4, 5e-5, 2e-5], save='yes')
#plt_tvstau('rad_infall2e-3_RK4_corr', 2e-3, 'no')
#plt_tvstau('rad_infall1e-3_RK4_corr', 1e-3, 'no')
#plt_tvstau('rad_infall5e-4_RK4_corr', 5e-4, 'no')
#plot_orbit('rad_infall1e-3_RK4','Radial Infall', 'upper right', fold + 'radial_infall.eps')


''' Crazy infalls '''
#plot_orbit('infall','', 'upper right', fold + 'infall1.eps')
#plot_orbit('infall2','', 'upper right', fold + 'infall2.eps')
#plot_orbit('volevi', '', 'upper right', fold + 'volevi.eps')
#check_stability_infall2(hs, 'yes')
#check_time_infall2    ([2e-3, 1e-3, 4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5, 1e-5], 'no')
#check_stability_volevi([2e-3, 1e-3, 8e-4, 4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5, 1e-5], 'no')

''' Circular Orbits '''
#plot_orbit('circular_orbit3', '', 'upper left', '')
#check_circular('circular_orbit3', 1e-3, fold + 'circ_res.png')
#check_circular('circular_orbit3_RK4_corrected', 1e-3, '')
#check_circular('circular_orbit4', 1e-4)


''' Precession '''
#plot_orbit('precession1', '', 'upper left', '../latex/Figures/chapter2/prec1.eps')
#plot_orbit('precession2', '', 'upper left', '../latex/Figures/chapter2/prec2.eps')


plt.show()
