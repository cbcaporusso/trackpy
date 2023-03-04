import numpy as np
import subprocess

def compute_displacement(file_path: str, dspl_dt, time: int) -> None:
    """
    Compute the displacement of the particles in the system.
    
    Parameters
    ----------
    file_path : str
        The path to the configuration file.
    time : int
        The time at which to compute the displacement.

    Returns
    -------
    None

    """

    process_status = subprocess.run(["/gpfs/projects/ub35/demian/chiral/displacement.sh", file_path, str(time), str(dspl_dt)])
