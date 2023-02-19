import numpy as np
from glob import glob


class Conf():
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

    @property
    def basepath(self):
        return self.filepath()

    def filepath(self, sub=None, time=None):
        """
        Return the filepath of the configuration file.

        Returns
        -------
        filepath : str
            The filepath of the configuration file.

        """
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

        data = np.loadtxt(self.filepath('hexatic', time))

        if data.shape[1] < 6:
            raise ValueError("Hexatic data file has incorrect shape")

        pos = data[:, [1, 2]]
        psi6 = data[:, [4, 5]]

        return pos, psi6

    def load_configuration(self, time):
        """
        Load the configuration for a given time.

        Parameters
        ----------
        time : int
            The time at which to load the configuration.

        Returns
        -------
        conf_data : numpy.ndarray
            The configuration data.

        """

        conf_data = np.loadtxt(self.filepath('trj', time))
        
        if conf_data.shape[1] < 3:
            raise ValueError("Configuration data file has incorrect shape")

        pos = conf_data[:, [1, 2]]
        
        return pos