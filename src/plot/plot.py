# import 

import sys, os, subprocess
sys.path.append('/gpfs/projects/ub35/demian/chiral')

import numpy as np
import matplotlib.pyplot as plt

from src.traj.box import Box
from matplotlib.collections import EllipseCollection
from typing import Tuple
from glob import glob

LJ_TIMESTEP = 0.01

# TODO: use setuptools to install the package and import it


# TODO add all the mode and hexatic part
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

    rho = f"{rho:.3f}"

    file_path=f"N_{N}/sigma_5.0/omega_{omega}/T_{temp}/rho_{rho}_1"

    if time == 'last':
        time = max([int(f.split('.')[-1]) for f in glob(f"{file_path}/Trj/xyz.dump.*")])   
    elif time < 0:
        time_array = sorted([int(f.split('.')[-1]) for f in glob(f"{file_path}/Trj/xyz.dump.*")])[time:]
        for t in time_array:
            iter_fig, iter_ax = plot_configuration(N, temp, omega, rho, t, save, figsize, print_params)
            if t == time_array[-1]:
                return iter_fig, iter_ax
            plt.close(iter_fig)

    print("Elaborating image at time: ", time)

    if not isinstance(time, int) and time != 'last':
        raise ValueError("time must be an integer or 'last'")

    sim_box = Box() 
    sim_box.read_box_size_from_file(f"{file_path}/Trj/xyz.dump.{time}")

    if not os.path.exists(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic"):
        print("No hexatic data found for the specified time, computing now...")
        compute_local_hexatic(file_path, time)
        
    data = np.loadtxt(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic")
    
    pos         = data[:,[1,2]]
    psi6_re     = data[:,4]
    psi6_im     = data[:,5]
    hex_args    = np.arctan2(psi6_im,psi6_re)

    # Plot the data.
    fig, ax = plt.subplots(figsize=figsize)
    set_lim(ax, sim_box)
    ec = add_coloured_collection(pos, hex_args, ax)
    
    if print_params:
        plot_params(ax, N, temp, omega, rho, time)

    # check if save is True or a string
    if save:
        if isinstance(save, str):
            fig.save(save)
        else:
            out_dir = f"snapshots/N_{N}/sigma_5.0/omega_{omega}/T_{temp}/rho_{rho}/hexatic"
            os.makedirs(out_dir, exist_ok=True)
            fig.savefig(f"{out_dir}/xyz.dump.{time}.png",
                        dpi=300, bbox_inches='tight')

    return fig, ax


def plot_configuration_from_path(file_path, savedir, time='last',  
                    mode='hexatic', figsize=(2.0,2.0), print_params=False):
    """
    Plot the configuration of the system at a given time.

    Parameters
    ----------
    file_path : str
        Path to the simulation data.
    time : int, optional
        Time at which to plot the configuration. If 'last' (default),
        the last time step is used.
    save : str or bool, optional
        If False (default) do not save a png of the configuration.
        If True, save the plot to the default file location. If a string,
        save the plot to the specified file location.
    mode : str, optional
        The mode of the plot. Can be 'hexatic' (default) or 'density' TODO
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

    N, temp, omega, rho = extract_params_from_path(file_path)
    # TODO add a better file saving system
    outdir = f"{savedir}/N_{N}/sigma_5.0/omega_{omega}/T_{temp}/rho_{rho}"
    fig, ax = plot_configuration(N, temp, omega, rho, time, savedir, figsize, print_params)


def plot_slab_configuration(temp, omega, rho, time='last', 
                            hexatic=False, dspl=None,
                            save=False, print_params=False):
    """
    Plot the configuration of a slab of the system at a given time.
    
    Parameters
    ----------
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

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.  
    ax : matplotlib.axes.Axes
        The axes object.

    """

    rho = f"{rho:.3f}"

    file_path=f"sub_projects/slab/sintetic_slab_pressure/T_{temp}/omega_{omega}/rho_{rho}_1"

    if time == 'last':
        time = max([int(f.split('.')[-1]) for f in glob(f"{file_path}/Trj/xyz.dump.*")])   
    elif time < 0:
            time_array = sorted([int(f.split('.')[-1]) for f in glob(f"{file_path}/Trj/xyz.dump.*")])[time:]
            for t in time_array:
                iter_fig, iter_ax = plot_slab_configuration(temp, omega, rho, t, hexatic, dspl, save, print_params)
                if t == time_array[-1]:
                    return iter_fig, iter_ax
                plt.close(iter_fig)
        
    print("Elaborating image at time: ", time)
    
    sim_box = Box()
    sim_box.read_box_size_from_file(f"{file_path}/Trj/xyz.dump.{time}")

    size_ratio = sim_box.lx / sim_box.ly

    fig, ax = plt.subplots(figsize=(2.0,2.0/size_ratio))

    if hexatic:
        if dspl is not None:
            raise ValueError("Cannot plot hexatic and displacement at the same time")

        if not os.path.exists(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic"):
            print("No hexatic data found for the specified time, computing now...")
            compute_local_hexatic(f"{file_path}", time)
        data = np.loadtxt(f"{file_path}/local_hexatic/xyz.dump.{time}.hexatic")
        pos         = data[:,[1,2]]
        psi6_re     = data[:,4]
        psi6_im     = data[:,5]
        hex_args    = np.arctan2(psi6_im,psi6_re)
        ec = add_coloured_collection(pos, hex_args, ax)
    
    else:
        if not os.path.exists(f"{file_path}/Dspl_dt_{dspl}/xyz.dump.{time}.dspl"):
            raise FileNotFoundError("No dspl data found for the specified time")

        data = np.loadtxt(f"{file_path}/Dspl_dt_{dspl}/xyz.dump.{time}.dspl")
        pos = data[:,[1,2]]
        dspl_mod = np.sqrt(data[:,3]**2 + data[:,4]**2)
        ec = add_coloured_collection(pos, dspl_mod, ax)

    set_lim(ax, sim_box)

    if print_params:
        plot_params(ax, '', temp, omega, rho, time)

    if save:
        if isinstance(save, str):
            fig.save(save)
        else:
            if hexatic:
                out_dir = f"snapshots/slab/T_{temp}/omega_{omega}/rho_{rho}/hexatic"
            else:
                out_dir = f"snapshots/slab/T_{temp}/omega_{omega}/rho_{rho}/displ"
            os.makedirs(out_dir, exist_ok=True)
            fig.savefig(f"{out_dir}/xyz.dump.{time}.png", dpi=300, bbox_inches='tight')

    return fig, ax 

def add_coloured_collection(pos, cmap, ax):
    """
    Add a hexatic collection to the axes.

    Parameters
    ----------
    pos : numpy.ndarray
        The positions of the particles.
    cmap : numpy.ndarray
        The array which will be converted to a colour map.
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
    ec.set_array(cmap)
    ax.add_collection(ec)
    
    return ec

def set_lim(ax, box: Box, shrink=0.0):
    """
    Set the limits of the axes.

    Parameters
    ----------
    box : box.Box
        The box object.
    ax : matplotlib.axes.Axes
        The axes object.
    shrink : float, optional
        The amount to shrink the axes by. Default is 0.05.

    Returns
    -------
    None

    """

    ax.set_xlim(0.0 + shrink * box.lx, box.lx - shrink * box.lx)
    ax.set_ylim(0.0 + shrink * box.ly, box.ly - shrink * box.ly)


def compute_local_hexatic(file_path, time):
    """
    Compute the local hexatic order parameter for a given time, calling an external bash script.

    Parameters
    ----------
    N : int
        Number of particles in the system.
    temp : float
        Temperature of the system.
    omega : float
        Frequency of the chiral force.
    rho : float 
        Density of the system.
    time : int  
        Time at which to compute the local hexatic order parameter.

    Returns
    -------
    None

    """
    process_status = subprocess.run(["/gpfs/projects/ub35/demian/chiral/elaborate.sh", file_path, str(time)])
    return process_status

def compute_coarsegrained_displacement():
    raise NotImplementedError

def plot_params(ax, N, temp, omega, rho, time):
    ax.text(0.05, 0.92, 
        f"$N = {N}, \\rho = {rho}, \\omega = {omega}, T = {temp}$",
        transform=ax.transAxes, fontsize=5)
    ax.text(0.05, 0.08,
        f"$t = {int(time * LJ_TIMESTEP)}$",
        transform=ax.transAxes, fontsize=5)

def extract_params_from_path(filepath: str) -> Tuple[int, float, float, float]:
    """
    Extract the parameters from a file path.

    Parameters
    ----------
    filepath : str
        The path to the file.

    Returns
    -------
    N : int
        Number of particles in the system.
    temp : float
        Temperature of the system.
    omega : float
        Frequency of the chiral force.
    rho : float
        Density of the system.
    time : int
        Time at which to compute the local hexatic order parameter.

    """

    N = int(filepath.split('/')[find_param_index_in_path(filepath, 'N_')].split('_')[1])
    temp = float(filepath.split('/')[find_param_index_in_path(filepath, 'T_')].split('_')[1])
    omega = float(filepath.split('/')[find_param_index_in_path(filepath, 'omega_')].split('_')[1])
    rho = float(filepath.split('/')[find_param_index_in_path(filepath, 'rho_')].split('_')[1])

    return N, temp, omega, rho


def find_param_index_in_path(filepath: str, param: str) -> int:
    """
    Find the index of a parameter in a file path.

    Parameters
    ----------
    filepath : str
        The path to the file.
    param : str
        The parameter to find.

    Returns
    -------
    index : int
        The index of the parameter in the file path.

    """

    list = filepath.split('/')
    # check partial matching pattern 'rho_*' in list
    # and extract the index of the first match
    try:
        index = next(i for i, s in enumerate(list) if param in s)
    except StopIteration:
        print("No rho parameter found in path.")
        exit()

    return index
