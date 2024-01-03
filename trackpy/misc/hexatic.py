import numpy as np
import subprocess
import os


def get_main_dir():
    package_path = os.path.dirname(__file__)
    package_path = os.path.abspath(os.path.join(package_path, ".."))
    package_path = os.path.abspath(os.path.join(package_path, ".."))
    return package_path


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
    # TODO: use subrocess to call the script in the src/bin subfolder
    script_path = os.path.join(get_main_dir(), "ext_scripts")
    print(script_path)
    process_status = subprocess.run(
        [f"{script_path}/elaborate.sh", script_path, file_path, str(time)]
    )
    return process_status


def global_hex_parameter(local_psi6: np.ndarray) -> np.ndarray:
    """
    Compute the global hexatic order parameter.

    Parameters
    ----------
    local_psi6 : numpy.ndarray
        The local hexatic order parameter.
        It is an array of shape (N, 2), where N is the number of particles,
        and the second dimension contains the real and imaginary parts
        of the local hexatic order parameter.

    Returns
    -------
    psi6 : float
        The global hexatic order parameter.

    """

    psi6_re = np.mean(local_psi6[:, 0])
    psi6_im = np.mean(local_psi6[:, 1])

    return np.array([psi6_re, psi6_im])


def global_hex_modulus(local_psi6: np.ndarray) -> float:
    """
    Compute the global hexatic order parameter modulus.

    Parameters
    ----------
    local_psi6 : numpy.ndarray
        The local hexatic order parameter.
        It is an array of shape (N, 2), where N is the number of particles,
        and the second dimension contains the real and imaginary parts
        of the local hexatic order parameter.

    Returns
    -------
    psi6 : float
        The global hexatic order parameter modulus.

    """

    global_psi6 = global_hex_parameter(local_psi6)
    return np.sqrt(global_psi6[0] ** 2 + global_psi6[1] ** 2)
