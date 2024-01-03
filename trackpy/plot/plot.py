from glob import glob

import sys
import numpy as np

from matplotlib.collections import EllipseCollection
from ..proj.chiral import Conf
from ..traj.box import Box
#from lib.utils import compute_local_hexatic, extract_params_from_path

from ..statistic import cg_field_to_cells
from ..misc.hexatic import global_hex_parameter

sys.path.append('/gpfs/projects/ub35/demian/chiral')

LJ_TIMESTEP = 0.01

# TODO: use setuptools to install the package and import it
# TODO add all the mode and hexatic part

def plot_configuration(conf: Conf, time, ax):
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

    if not isinstance(time, int) or time < 0:
        raise ValueError("time must be a positive integer")
      
    pos = conf.load_configuration(time)
    
    sim_box = Box() 
    sim_box.read_box_size_from_file(conf.filepath(sub="trj", time=time)) 
    set_lim(ax, sim_box)
    
    ec = add_coloured_collection(pos, np.ones(pos.shape[0]), ax)
    
    return ec
     

def plot_hexatic_arg(conf: Conf, time, ax):

    if not isinstance(time, int) or time < 0:
        raise ValueError("time must be a positive integer")
    
    pos, psi6 = conf.load_hexatic(time)
    
    sim_box = Box() 
    sim_box.read_box_size_from_file(conf.filepath(sub="trj", time=time))
    set_lim(ax, sim_box)
    
    ec = add_coloured_collection(pos, np.tanh(psi6[1], psi6[0]), ax)
    ec.set_cmap('hsv')
    
    return ec


def plot_hexatic_mod(conf: Conf, time, ax):

    if not isinstance(time, int) or time < 0:
        raise ValueError("time must be a positive integer")
    
    pos, psi6 = conf.load_hexatic(time)
    
    sim_box = Box() 
    sim_box.read_box_size_from_file(conf.filepath(sub="trj", time=time))
    set_lim(ax, sim_box)
    
    ec = add_coloured_collection(pos, np.sqrt(psi6[0]**2 + psi6[1]**2), ax)
    ec.set_cmap('jet')
    
    return ec

def plot_hexatic_proj(conf: Conf, time, ax):

    if not isinstance(time, int) or time < 0:
        raise ValueError("time must be a positive integer")
    
    pos, psi6 = conf.load_hexatic(time)
    
    sim_box = Box() 
    sim_box.read_box_size_from_file(conf.filepath(sub="trj", time=time))
    set_lim(ax, sim_box)

    # here we compute the projection 
    psi6_glob = global_hex_parameter(psi6)
    psi6_proj = psi6_glob[0] * psi6[:,0] + psi6_glob[1] * psi6[:,1]
    psi6_proj /= np.sqrt(psi6_glob[0] ** 2 + psi6_glob[1] ** 2)
    psi6_proj /= np.sqrt(psi6[:,0] ** 2 + psi6[:,1] ** 2)
    
    ec = add_coloured_collection(pos, psi6_proj, ax)
    ec.set_cmap('jet')

    return ec

# def plot_configuration_from_path(file_path, savedir, time='last',  
#                     mode='hexatic', figsize=(2.0,2.0), print_params=False):
#     """
#     Plot the configuration of the system at a given time.

#     Parameters
#     ----------
#     file_path : str
#         Path to the simulation data.
#     time : int, optional
#         Time at which to plot the configuration. If 'last' (default),
#         the last time step is used.
#     save : str or bool, optional
#         If False (default) do not save a png of the configuration.
#         If True, save the plot to the default file location. If a string,
#         save the plot to the specified file location.
#     mode : str, optional
#         The mode of the plot. Can be 'hexatic' (default) or 'density' TODO
#     figsize : tuple, optional
#         The size of the figure in inches. Default is (2.0,2.0).
#     print_params : bool, optional
#         If True, print the parameters of the simulation in the plot. Default is False.

#     Returns
#     -------
#     fig : matplotlib.figure.Figure
#         The figure object.
#     ax : matplotlib.axes.Axes
#         The axes object.
#     """

#     N, temp, omega, rho = extract_params_from_path(file_path)
#     # TODO add a better file saving system
#     outdir = f"{savedir}/N_{N}/sigma_5.0/omega_{omega}/T_{temp}/rho_{rho}"
#     fig, ax = plot_configuration(N, temp, omega, rho, time, savedir, figsize, print_params)


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

    ec = EllipseCollection(
        1, 1, 0, units='xy', offsets=pos,
        transOffset=ax.transData)

    #ec.set_alpha(0.50)
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


def compute_coarsegrained_displacement():
    raise NotImplementedError


def plot_params(ax, N, temp, omega, rho, time):
    ax.text(0.05, 0.92, 
        f"$N = {N}, \\rho = {rho}, \\omega = {omega}, T = {temp}$",
        transform=ax.transAxes, fontsize=5)
    ax.text(0.05, 0.08,
        f"$t = {int(time * LJ_TIMESTEP)}$",
        transform=ax.transAxes, fontsize=5)


def add_configuration(pos, box, time, pmap, ax):
    """
    Add a configuration to the axes.

    Parameters
    ----------
    conf : numpy.ndarray
        The configuration to plot.
    ax : matplotlib.axes.Axes
        The axes object.

    Returns
    -------
    None

    """

    set_lim(ax, sim_box)
    ec = add_coloured_collection(pos, pmap, ax)

    return ec

def plot_droplet(ax, conf: Conf, time, dspl):
    """"""

    pos, displ = conf.load_displacement(time, dspl)
    
    # sim_box = Box.from_dump(conf.filepath(sub="trj", time=time))
    # set_lim(ax, sim_box)

    ec = add_coloured_collection(pos, np.sqrt(displ[:,0]**2 + displ[:,1]**2), ax)
    
    return ec

def plot_droplet_with_disp_field(ax, conf: Conf, time, dspl, cell_size = 2.5, scale=0.05):
    """"""

    pos, displ = conf.load_displacement(time, dspl)
    
    # sim_box = Box.from_dump(conf.filepath(sub="trj", time=time))
    # set_lim(ax, sim_box)

    cells, field_x = cg_field_to_cells(pos, cell_size, displ[:,0])
    cells, field_y = cg_field_to_cells(pos, cell_size, displ[:,1])

    ec = add_coloured_collection(pos, np.sqrt(displ[:,0]**2 + displ[:,1]**2), ax)
    qs = ax.quiver(cells[:,0], cells[:,1], field_x, field_y, scale=scale, color='black')

    return ec, qs 
    
