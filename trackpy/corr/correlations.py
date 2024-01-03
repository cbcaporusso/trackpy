import numpy as np
from .src.box import Box

def hex_correlations(psi_6, box):
    """
    Returns the hexatic correlations 

    Parameters
    ----------
    psi_6 : array
        The hexatic order parameter
    box : Box object
        
    Returns
    -------
    c : array
        The hexatic correlations
    """
    
    assert psi_6.ndim == 2, "psi_6 must be a 2D array"

    box_size = np.max(box.box_size)

    # compute the distance matrix between all particles
    distance_matrix = np.sqrt(np.sum((box.points[:, None, :] - box.points[None, :, :])**2, axis = -1))
    distance_matrix = np.minimum(distance_matrix, box_size - distance_matrix)