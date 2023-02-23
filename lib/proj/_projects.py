import os
from abc import ABC, abstractmethod
from glob import glob

import numpy as np

from misc.hexatic import compute_local_hexatic


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

        pos = conf_data[:, [1, 2]]
        
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

