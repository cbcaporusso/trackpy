import numpy as np
import logging
from ..traj.box import Box

""" 
This module contains functions for making initial configurations for simulations.
"""

def place_hexagonally(num_part: int, box: Box):
    """
    Place particles in an hexagonal lattice.

    Parameters
    ----------
    num_part : int
        The number of particles to place.
    box: Box
        The simulation box.
    
    Returns
    -------
    positions : np.ndarray
        The positions of the particles.
    """

    if not isinstance(num_part, int):
        raise TypeError("The number of particles must be an integer.")
    
    if num_part < 0:
        raise ValueError("The number of particles must be positive.")
    
    if not isinstance(box, Box):
        raise TypeError("The simulation box must be an instance of Box.")
        
    # Define the basis vectors
    a = 1.0
    a1 = np.array([a, 0])
    a2 = np.array([a/2, a*np.sqrt(3)/2])

    # Determine the number of rows and columns of the lattice
    nrows = int(np.sqrt(num_part / (4*np.sqrt(3)/3)))
    ncols = int(np.ceil(num_part / nrows))

    # Generate the positions of the particles
    positions = []
    for i in range(-nrows, nrows):
        for j in range(-ncols, ncols):
            x = i*a1[0] + j*a2[0]
            y = i*a1[1] + j*a2[1]
            # Check if the particle is within the simulation box
            if box.isinbox([x, y]):
                positions.append([x, y, 0.0])

    # Return the positions of the particles
    return np.array(positions)


def remove_particles(positions, box: Box, r_cut):
    
    # remove the particles outside the box
    
    # Compute the distance of each particle from the center
    r0 = box.center
    distances = box.distance_pbc(positions, r0)

    # Select the particles that are closer than r_cut
    mask = distances < r_cut

    # Return the positions of the particles
    return positions[mask]


def add_gas(pos, box: Box, target_density, log=False):
    # add gas particles in random pos in the box
    # checking that they do not overlap with the particles in pos
    MAX_PARTICLES = 1024 ** 2
    gas = np.zeros((MAX_PARTICLES, 3))

    density = compute_density(pos, box)
    
    count = 0
    while density < target_density:
        if log:
            logging.info("density: {:.3f}".format(density) +
                         " target: {:.3f}".format(target_density))
        # here we define boundary
        # such as there is a 0.5 distance
        # between the edges of the box and the particles
        boundary = 0.5
        
        # this generates a random position
        # between [0, L - 2 * sigma] and then 
        # shifts it by 0.5 * sigma  
        new_pos = \
            np.random.rand(3)*[box.lx - 2 * boundary,
                               box.ly - 2 * boundary,
                               0.0]
        new_pos += [boundary, boundary, 0.0]

        if not overlap(pos, new_pos) and \
                not overlap(gas[:count, :], new_pos):
            gas[count, :] = new_pos
            density += np.pi/4/box.volume
            count += 1
    return gas[:count, :]


def overlap(pos, new_pos):
    # check if the new particle overlaps with the existing ones
    for pp in pos:
        if np.linalg.norm(pp-new_pos) <= 1.5:
            return True
    return False


def compute_density(pos, box: Box):
    # compute the density of the system
    print(pos)
    N = len(pos[:, 0])
    return N*np.pi/4/box.volume


def count_particles(pos):
    return len(pos[:, 0])

