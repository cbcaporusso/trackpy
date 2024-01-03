from typing import Tuple
import numpy as np
from ..proj._projects import Subproject
import pandas as pd


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

    N = int(
        filepath.split("/")[_find_param_index_in_path(filepath, "N_")].split("_")[1]
    )
    temp = float(
        filepath.split("/")[_find_param_index_in_path(filepath, "T_")].split("_")[1]
    )
    omega = float(
        filepath.split("/")[_find_param_index_in_path(filepath, "omega_")].split("_")[1]
    )
    rho = float(
        filepath.split("/")[_find_param_index_in_path(filepath, "rho_")].split("_")[1]
    )

    return N, temp, omega, rho


def _find_param_index_in_path(filepath: str, param: str) -> int:
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

    split_path_list = filepath.split("/")
    # check partial matching pattern 'rho_*' in list
    # and extract the index of the first match
    try:
        index = next(i for i, s in enumerate(split_path_list) if param in s)
    except StopIteration:
        print("No parameter found in path.")
        exit()

    return index


def distance_pbc(x0, x1, L) -> np.ndarray:
    """Computes the distances using periodic boundary
    conditions (PBC)"""
    delta = x0 - x1
    delta = np.where(delta > 0.5 * L, delta - L, delta)
    delta = np.where(delta < -0.5 * L, delta + L, delta)
    return delta


def norm(x: np.ndarray) -> float:
    return np.sqrt((x**2).sum(axis=-1))


def retrive_particle_trj_from_id(
    conf: Subproject, atomid: int, time_start=None, time_end=None
):
    """
    Extract the trajectory of a particle from its id.

    Parameters
    ----------
    conf : Subproject
        The configuration object.
    atomid : int
        The id of the particle.
    time_start : int, optional
        The starting time of the trajectory.
    time_end : int, optional
        The ending time of the trajectory.

    Returns
    -------
    trj : numpy.ndarray
        The trajectory of the particle.
    """

    time_array = np.array(conf.find_configuration_times("*"))

    if time_start is not None:
        time_array = time_array[time_array >= time_start]
    if time_end is not None:
        time_array = time_array[time_array <= time_end]

    trj = np.zeros((len(time_array), 2))
    for i, time in enumerate(time_array):
        file_path = conf.filepath("trj", time)
        df = pd.read_csv(file_path, skiprows=8, sep=" ")
        df = df[df["ITEM:"] == atomid]
        trj[i] = df[["ATOMS", "id"]].to_numpy()

    return trj
