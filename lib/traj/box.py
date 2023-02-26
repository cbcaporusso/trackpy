import numpy as np 
from utils.lammps import lammps_boxsize_parser

# In future this class should be changed.
# The box size should be a 4 element array,
# representing the simulation box
# TODO: 
# - change initializer to account for the 4 element array

class Box:
    """
    Class to store the size of the simulation box


    Methods
    -------
    read_box_size_from_file(filename)
        read the size of the simulation box from a file given as a parameter


    """
    
    def __init__(self, lx, ly) -> None:
        self.lx = lx 
        self.ly = ly  

        self.L = np.array([lx, ly])
        

    def __str__(self) -> str:
        return f"lx = {self.lx}, ly = {self.ly}"


    # TODO: we could also add a method to read the box size from a a subproject object
    @staticmethod
    def from_dump(filename: str) -> None:
        """
        Read the size of the simulation box from a file given as a parameter.

        Parameters
        ----------
        filename : str
            The name of the file to read the box size from.

        Returns
        -------
        box : Box
            The box object.

        """

        lx, ly = lammps_boxsize_parser(filename)
        return Box(lx, ly)


    def distance_pbc(self, x0, x1) -> np.ndarray:
        """
        Compute the distance between two points in a periodic box.

        Parameters
        ----------
        x0 : np.ndarray
            The first point.
        x1 : np.ndarray
            The second point.

        Returns
        -------
        delta : np.ndarray
            The distance between the two points.

        """
        
        delta = x0 - x1
        delta = np.where(delta > 0.5 * self.L, delta - self.L, delta)
        delta = np.where(delta < - 0.5 * self.L, delta + self.L, delta)
        return delta


    @staticmethod
    def norm(x: np.ndarray) -> float:
        """
        Compute the norm of a vector.
        
        Parameters
        ----------
        x : np.ndarray
            The vector.
        
        Returns
        -------
        norm : float
            The norm of the vector.
        
        """

        return np.sqrt((x ** 2).sum(axis=-1))
