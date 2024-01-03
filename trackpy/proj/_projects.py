import os
from abc import ABC, abstractmethod
from glob import glob

import numpy as np

from ..misc.hexatic import compute_local_hexatic
from ..misc.displacement import compute_displacement
from ..utils.lammps import lammps_header_parser

# TODO - add an autmated way to find the indexes


class Subproject(ABC):
    """
    Abstract class for subprojects.
    """

    @abstractmethod
    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def filepath(self, sub=None, time=None):
        """
        Return the filepath of the configuration file.

        Returns
        -------
        filepath : str
            The filepath of the configuration file.

        """
        pass

    @property
    def basepath(self):
        return self.filepath()

    def find_configuration_times(self, last_vals):
        """
        Find the configuration files for the system.

        Returns
        -------
        confs : list
            A list of the configuration files.

        """

        if isinstance(last_vals, str) and last_vals == "*":
            return sorted([int(f.split('.')[-1])
                           for f in glob(self.filepath('trj', '*'))])

        if isinstance(last_vals, int):
            if last_vals > 0:
                raise ValueError("lastvals must be negative or *")
            return sorted([int(f.split('.')[-1])
                           for f in glob(self.filepath('trj', '*'))])[last_vals:]   

    def load_configuration(self, time):
        """
        Load the configuration at a given time.

        Parameters
        ----------
        time : int
            The time at which to load the configuration.

        Returns
        -------
        conf : numpy.ndarray
            The configuration.

        """

        conf_data = np.loadtxt(self.filepath('trj', time), skiprows=9)
        
        if conf_data.shape[1] < 3:
            raise ValueError("Configuration data file has incorrect shape")

        pos_index = [lammps_header_parser(self.filepath('trj', time), "x"),
                     lammps_header_parser(self.filepath('trj', time), "y")]
        
        pos = conf_data[:, pos_index]
        
        return pos

    def load_hexatic(self, time):
        """
        Load the hexatic data for a given time.

        Parameters
        ----------
        time : int
            The time at which to load the hexatic data.

        Returns
        -------
        hex_data : numpy.ndarray
            The hexatic data.

        """

        if not os.path.exists(self.filepath(sub="hexatic", time=time)):
            print("No hexatic data found for the specified time, computing now...")
            compute_local_hexatic(self.basepath, time)

        data = np.loadtxt(self.filepath('hexatic', time))

        if data.shape[1] < 6:
            raise ValueError("Hexatic data file has incorrect shape")

        pos = data[:, [1, 2]]
        psi6 = data[:, [4, 5]]

        return pos, psi6

    def load_displacement(self, time, dspl_dt):
        """
        Load the displacement data for a given time.

        Parameters
        ----------
        time : int
            The time at which to load the displacement data.
        dspl_dt : int
            The time interval between configurations.

        Returns
        -------
        dspl_data : numpy.ndarray
            The displacement data.

        """

        if not os.path.exists(self.filepath(sub=f"dspl_{dspl_dt}", time=time)):
            print("No displacement data found for the specified time, computing now...")
            compute_displacement(self.basepath, dspl_dt, time)

        data = np.loadtxt(self.filepath(f"dspl_{dspl_dt}", time), skiprows=9)

        if data.shape[1] < 4:
            raise ValueError("Displacement data file has incorrect shape")

        pos_index = [lammps_header_parser(self.filepath(f'dspl_{dspl_dt}', time), "x"),
                     lammps_header_parser(self.filepath(f'dspl_{dspl_dt}', time), "y")]
        vel_index = [lammps_header_parser(self.filepath(f'dspl_{dspl_dt}', time), "vx"),
                     lammps_header_parser(self.filepath(f'dspl_{dspl_dt}', time), "vy")]

        pos = data[:, pos_index]
        dspl = data[:, vel_index]

        return pos, dspl

    @staticmethod
    def folder_structure(sub: str = None, time : int = None):
        """
        Return the folder structure for the subproject.
        This is the standard folder structure, if subclasses
        have a different folder structure, this method should
        be overwritten.

        Parameters
        ----------
        sub : str
            The subfolder to return the folder structure for.
        time : int, optional
            The time of the configuration.

        Returns
        -------
        out_dir : str
            The folder structure for the subproject.

        """

        if sub == "trj":
            out_dir = "/Trj/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}"
            else:
                return out_dir

        elif sub == "hexatic":
            out_dir = "/local_hexatic/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}.hexatic"
            else:
                return out_dir

        elif sub == None:
            return ''

        elif isinstance(sub, str) and "dspl_" in sub:
            dt_dspl = sub.split("_")[1]
            out_dir = f"/Dspl_dt_{dt_dspl}/"
            if time is not None:
                return out_dir + f"xyz.dump.{time}.dspl"
            else:
                return out_dir

        else:
            raise ValueError("Unknown subproject")    
