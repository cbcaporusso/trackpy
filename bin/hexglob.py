import numpy as np
import os

from lib.proj.chiral import Conf
from lib.misc.hexatic import (
    compute_local_hexatic, global_hex_modulus)

if __name__ == "__main__":

    npart = 262144
    temp = 0.35
    omegas = [1.0]
    rho = 0.700

    for omega in omegas:
        conf = Conf(npart, temp, omega, rho)
        glob_array = []
        for time in conf.find_configuration_times(-10):
            # check if the file exists
            if not os.path.isfile(conf.filepath(sub="hexatic", time=time)):
                print("No hexatic data found for the specified time, computing now...")
                compute_local_hexatic(conf.filepath(), time)

            pos, hex_data = conf.load_hexatic(time)
            glob_psi6 = global_hex_modulus(hex_data)
            glob_array.append(glob_psi6)

        glob_psi6_avg = np.mean(glob_array)
