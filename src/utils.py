import numpy as np
import subprocess

from typing import Tuple

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

    N = int(filepath.split('/')[_find_param_index_in_path(filepath, 'N_')].split('_')[1])
    temp = float(filepath.split('/')[_find_param_index_in_path(filepath, 'T_')].split('_')[1])
    omega = float(filepath.split('/')[_find_param_index_in_path(filepath, 'omega_')].split('_')[1])
    rho = float(filepath.split('/')[_find_param_index_in_path(filepath, 'rho_')].split('_')[1])

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

    list = filepath.split('/')
    # check partial matching pattern 'rho_*' in list
    # and extract the index of the first match
    try:
        index = next(i for i, s in enumerate(list) if param in s)
    except StopIteration:
        print("No rho parameter found in path.")
        exit()

    return index