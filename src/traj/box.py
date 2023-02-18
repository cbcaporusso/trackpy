import numpy as np 

class Box:
    """
    Class to store the size of the simulation box


    Methods
    -------
    read_box_size_from_file(filename)
        read the size of the simulation box from a file given as a parameter


    """
    
    def __init__(self) -> None:
        lx : float
        ly : float 
        L : np.array

    def __str__(self) -> str:
        return f"{self.L}"

    def read_box_size_from_file(self, filename: str) -> None:
        """
        read the size of the simulation box from a file given as a parameter 

        """
        with open(filename) as fp:
            for i, line in enumerate(fp):
                if i==5:
                    self.lx = np.fromstring(line, sep=' ')[1]
                if i==6:
                    self.ly = np.fromstring(line, sep=' ')[1]
                    self.L = np.array([self.lx,self.ly])
                    break

    def update(self):
        self.L = np.array([self.lx,self.ly])
