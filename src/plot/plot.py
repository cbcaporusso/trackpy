import sys, os
sys.path.append('/gpfs/projects/ub35/demian/chiral')

import numpy as np
import matplotlib.pyplot as plt
import src.box as box
import subprocess
from matplotlib.collections import EllipseCollection
from glob import glob

# TODO: use setuptools to install the package and import it

def plot_configuration(N, temp, omega, rho, time='last', save=False, 
                    figsize=(2.0,2.0), print_params=False):
    """
    Plot the configuration of the system at a given time.

    Parameters
    ----------
    N : int
        Number of particles.
    temp : float
        Temperature of the system.
    omega : float
        Frequency of the chiral force.
    rho : float
        Density of the system.
    time : int, optional
        Time at which to plot the configuration. If 'last' (default),
        the last time step is used.
    save : str or bool, optional
        If False (default) do not save a png of the configuration. 
        If True, save the plot to the default file location. If a string,
        save the plot to the specified file location.
    figsize : tuple, optional
        The size of the figure in inches. Default is (2.0,2.0).
    print_params : bool, optional
        If True, print the parameters of the simulation in the plot. Default is False.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    ax : matplotlib.axes.Axes
        The axes object.
    """
    # Load the data.

    file_path=f"N_{N}/sigma_5.0/omega_{omega}/T_{temp}/rho_{rho:.3f}_1"
    print(file_path)

    if time == 'last':
        #print(glob(f"{file_path}/Trj/xyz.dump.*"))
        time = max([int(f.split('.')[-1]) for f in glob(f"{file_path}/Trj/xyz.dump.*")])   
        print(time)
        #exit()

    if not isinstance(time, int) and time != 'last':
        raise ValueError("time must be an integer or 'last'")

    sim_box = box.Box() 
    sim_box.read_box_size_from_file(f"{file_path}/Trj/xyz.dump.{time}")

    if not os.path.exists(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic"):
        print("No hexatic data found for the specified time, computing now...")
        subprocess.run(["bash", "/gpfs/projects/ub35/demian/chiral/elaborate.bash", f"{N}", f"{temp}", f"{omega}", f"{rho:.3f}", f"{time}"])
        
    data = np.loadtxt(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic")
    
    pos         = data[:,[1,2]]
    psi6_re     = data[:,4]
    psi6_im     = data[:,5]
    hex_args    = np.arctan2(psi6_im,psi6_re)

    # Plot the data.
    fig, ax = plt.subplots(figsize=figsize)
    ec = add_hexatic_collection(pos, hex_args, ax)

    shrink_factor = 0.5

    ax.set_xlim(0, sim_box.lx * shrink_factor) 
    ax.set_ylim(0, sim_box.ly * shrink_factor)
    
    if print_params:
        ax.text(0.05, 0.95, f"$N = {N}, \\rho = {rho}, \\omega = {omega}, T = {temp}$", transform=ax.transAxes)

    # check if save is True or a string
    if save:
        if isinstance(save, str):
            fig.save(save)
        else:
            out_dir = f"snapshots/N_{N}/sigma_5.0/omega_{omega}/T_{temp}/rho_{rho:.3f}"
            os.makedirs(out_dir, exist_ok=True)
            fig.savefig(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic.png", dpi=300, bbox_inches='tight')

    return fig, ax

def plot_slab_configuration(temp, omega, rho):
    """
    Plot the configuration of a slab of the system at a given time.
    """

    file_path=f"sub_projects/slab/sintetic_slab_pressure/T_{temp}/omega_{omega}/rho_{rho:.3f}_1"

    size_ratio = sim_box.lx / sim_box.ly

    return fig, ax 

def add_hexatic_collection(pos, hex_args, ax):
    """
    Add a hexatic collection to the axes.

    Parameters
    ----------
    pos : numpy.ndarray
        The positions of the particles.
    hex_args : numpy.ndarray
        The hexatic arguments of the particles.
    ax : matplotlib.axes.Axes
        The axes object.

    Returns
    -------
    ec : matplotlib.collections.EllipseCollection
        The hexatic collection.

    """

    ax.set_xticks([])
    ax.set_yticks([])

    ec = EllipseCollection(1, 1, 0, 
    units='xy',
    offsets=pos,
    transOffset=ax.transData
    )

    ec.set_alpha(0.50)
    ec.set_cmap('hsv')
    ec.set_array(hex_args)
    ax.add_collection(ec)
    
    return ec

if __name__ == "__main__":
    plot_configuration(262144, 0.35, 3.0, 0.700, time='last', save=True, figsize=(2.0,2.0), print_params=True)