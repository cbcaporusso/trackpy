import numpy as np
import matplotlib.pyplot as plt

from output import LAMMPS_initialconf

import logging


def place_hexagonally(num_part, Lx, Ly):
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
            # Compute the position of the particle based on the row and column indices
            x = i*a1[0] + j*a2[0]
            y = i*a1[1] + j*a2[1]
            # Check if the particle is within the simulation box
            if x < Lx and y < Ly and x > 0 and y > 0:
                positions.append([x, y, 0.0])

    # Return the positions of the particles
    return np.array(positions)


def remove_particles(positions, Lx, Ly, r_cut):
    
    # remove the particles outside the box
    
    # Compute the distance of each particle from the center
    r0 = np.array([0.5*Lx, 0.5*Ly])
    distances = np.sqrt((positions[:,0]-r0[0])**2 + (positions[:,1]-r0[1])**2)

    # Select the particles that are closer than r_cut
    mask = distances < r_cut

    # Return the positions of the particles
    return positions[mask]


def print_LAMMPS_configurations(filename, positions, Lx, Ly):
    """print in output the configuration."""

    N = int(len(positions[:, 0]))
    
    with open(filename, 'w') as f:

        f.write("LAMMPS data file\n\n")
        f.write("{:d} atoms\n".format(N))
        # f.write("{:d} bonds\n\n".format(N))
        f.write("1 atom types\n")
        # f.write("1 bond types\n\n")

        f.write("0 {:f} xlo xhi\n".format(Lx))
        f.write("0 {:f} ylo yhi\n".format(Ly))
        f.write("-0.5 0.5 zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write("1 1\n\n")
        f.write("Atoms\n\n")

        for i in range(N):
            print(i)
            out = "{:d} 1 1 {:f} {:f} ".format(
                i+1, positions[i, 0], positions[i, 1])
            out += "0.0000 0 0 0\n"
            f.write(out)


def add_gas(pos, Lx, Ly, target_density, log=False):
    # add gas particles in random pos in the box
    # checking that they do not overlap with the particles in pos
    gas = np.zeros((10000, 3))
    density = compute_density(pos, Lx, Ly)
    count = 0
    while density < target_density:
        if log:
            logging.info("density: {:.3f}".format(density) +
                         " target: {:.3f}".format(target_density))
        # here we define boundary
        # such as there is a 0.5 distance
        # between the edges of the box and the particles
        boundary = 0.5
        new_pos = \
            np.random.rand(3)*[Lx - 2 * boundary, Ly - 2 * boundary, 0.0]
        new_pos += [boundary, boundary, 0.0]

        if not overlap(pos, new_pos) and \
                not overlap(gas[:count, :], new_pos):
            gas[count, :] = new_pos
            density += np.pi/4/Lx/Ly
            count += 1
    return gas[:count, :]


def compute_density(pos, Lx, Ly):
    # compute the density of the system
    print(pos)
    N = len(pos[:, 0])
    return N*np.pi/4/Lx/Ly


def count_particles(pos):
    return len(pos[:, 0])


def overlap(pos, new_pos):
    # check if the new particle overlaps with the existing ones
    for pp in pos:
        if np.linalg.norm(pp-new_pos) <= 1.5:
            return True
    return False


if __name__ == "__main__":
    # Define the number of particles
    num_part = 100000

    logging.basicConfig(level=logging.INFO)

    # Define the size of the simulation box
    Lx = 100
    Ly = 100

    initial_conf = LAMMPS_initialconf([Lx, Ly])
    
    # Define the radius of the circle
    r_cut = 20.0
    # Place the particles in an hexagonal lattice
    positions = place_hexagonally(num_part, Lx, Ly)
    # Remove the particles that are too far from the center
    cluster_pos = remove_particles(positions, Lx, Ly, r_cut)
    n_cluster = count_particles(cluster_pos)

    initial_conf.add_particles(cluster_pos, 1)

    # add gas particles
    gas_pos = add_gas(cluster_pos, Lx, Ly, 0.20, log=True)
    initial_conf.add_particles(gas_pos, 2)

    #tot_part = initial_conf.n_particles

    initial_conf.write(f"hex_with_gas/c_{n_cluster}.lmp")