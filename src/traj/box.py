import numpy as np 

class Box:

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

def distance_pbc(x0, x1, L)->np.ndarray:
    ''' Computes the distances using periodic boundary
    conditions (PBC)'''
    delta = x0 - x1
    delta = np.where(delta > 0.5 * L, delta - L, delta)
    delta = np.where(delta < - 0.5 * L, delta + L, delta)
    return delta

def norm(x: np.ndarray)->float:
    return np.sqrt((x ** 2).sum(axis=-1))
