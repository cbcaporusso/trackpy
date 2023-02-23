import os
import numpy as np
from glob import glob
from misc.hexatic import compute_local_hexatic
from lib.proj._projects import Subproject

# TODO: we could include the box class here


class Conf(Subproject):
    """
    A class to represent a configuration of particles.

    """

    def __init__(self, npart, temp, omega, rho):
        """
        Initialise a configuration object.

        Parameters
        ----------
        N : int
            The number of particles.
        temp : float
            The temperature of the system.
        omega : float
            The frequency of the system.
        rho : float
            The density of the system.

        Returns
        -------
        None

        """

        self.npart = npart
        self.temp = temp

        omega_decimals = str(round(omega - int(omega), 2))
        if omega_decimals == '0.0':
            self.omega = f"{omega:.1f}"
        else:
            self.omega = f"{omega:.2f}"

        self.rho = f"{rho:.3f}"

    def filepath(self, sub=None, time=None):
        base_dir = \
            f"N_{self.npart}/sigma_5.0/omega_{self.omega}/T_{self.temp}/rho_{self.rho}_1"

        if sub == "trj":
            out_dir = base_dir + "/Trj/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}"

        if sub == "hexatic":
            out_dir = base_dir + "/local_hexatic/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}.hexatic"

        return f"{base_dir}/"

class Disordered(Subproject):
    """
    A class to represent a configuration of particles.

    """

    def __init__(self, npart, temp, omega, rho):
        """
        Initialise a configuration object.

        Parameters
        ----------
        N : int
            The number of particles.
        temp : float
            The temperature of the system.
        omega : float
            The frequency of the system.
        rho : float
            The density of the system.

        Returns
        -------
        None

        """

        self.npart = npart
        self.temp = temp

        omega_decimals = str(round(omega - int(omega), 2))
        if omega_decimals == '0.0':
            self.omega = f"{omega:.1f}"
        else:
            self.omega = f"{omega:.2f}"

        self.rho = f"{rho:.3f}"

    def filepath(self, sub=None, time=None):
        base_dir = (f"N_{self.npart}/from_disordered/sigma_5.0/"
                    f"omega_{self.omega}/T_{self.temp}/rho_{self.rho}_1")

        if sub == "trj":
            out_dir = base_dir + "/Trj/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}"

        if sub == "hexatic":
            out_dir = base_dir + "/local_hexatic/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}.hexatic"

        return f"{base_dir}/"


class Slab(Subproject):
    def __init__(self, temp, omega, rho):
        """
        Initialise a slab object.

        Parameters
        ----------
        temp : float
            The temperature of the system.
        omega : float
            The frequency of the system.
        rho : float
            The density of the system.

        Returns
        -------
        None

        """

        self.temp = temp

        omega_decimals = str(round(omega - int(omega), 2))
        if omega_decimals == '0.0':
            self.omega = f"{omega:.1f}"
        else:
            self.omega = f"{omega:.2f}"

        self.rho = f"{rho:.3f}"

    def filepath(self, sub=None, time=None):
        base_dir = (f"sub_projects/slab/sintetic_slab_pressure/"
                    f"T_{self.temp}/omega_{self.omega}/rho_{self.rho}_1")
        if sub == "trj":
            out_dir = base_dir + "/Trj/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}"

        if sub == "hexatic":
            out_dir = base_dir + "/local_hexatic/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}.hexatic"

        return f"{base_dir}/"
