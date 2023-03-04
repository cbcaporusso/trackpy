import numpy as np


class LAMMPS_initialconf():
    def __init__(self, box: np.ndarray):
        self.lx = box[0]
        self.ly = box[1]

        self._num_types = 0
        self._num_particles = 0
        self._types_arr = []
        self._positions = []
        # self._open =

    def add_particles(self, positions: np.ndarray, part_type: int):
        if part_type not in self._types_arr:

            self._types_arr.append(part_type)
            self._num_types += 1

        self._positions.append(positions)
        self._num_particles += positions.shape[0]

    def _write_header(self, filename):
        with open(filename, 'w') as f:
            f.write("LAMMPS data file\n\n")
            f.write("{:d} atoms\n".format(self._num_particles))
            f.write("{:d} atom types\n\n".format(self._num_types))
            f.write("0 {:f} xlo xhi\n 0 {:f} ylo yhi\n -0.5 0.5 zlo zhi\n\n".format(self.lx, self.ly))
            f.write("Masses\n\n")
            for i, type in enumerate(self._types_arr):
                f.write("{:d} 1\n".format(i+1))
            f.write("\n")

    def write(self, filename: str):
        self._write_header(filename)
        with open(filename, 'a') as f:
            f.write("Atoms\n\n")
            for i, type in enumerate(self._types_arr):
                for j, pos in enumerate(self._positions[i]):
                    out = "{:d} {:d} 1 {:f} {:f} 0.0000 0 0 0\n".format(j+1, i+1, pos[0], pos[1])
                    f.write(out)

