import os
import numpy as np
from glob import glob
from ..misc.hexatic import compute_local_hexatic
from ._projects import Subproject

# TODO: we could include the box class here


class Conf(Subproject):
    """
    A class to represent a configuration of particles.

    """
    def __init__(self, npart, temp, omega, rho):
        pass

   
class Cluster(Subproject):
    raise NotImplementedError

