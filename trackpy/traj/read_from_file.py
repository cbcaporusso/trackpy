from typing import List
from interface import ReadInterface
import numpy as np
import os


class ReadFromFile(ReadInterface):
    ''' Implementation of the interface that read clusters from
    three files: Trj Radii and Labels'''

    POS_INDICES = (1,2)

    def __init__(self, path, extension=None) -> None:
        
        self.path = path
        self.extension = extension

    def read(self, time: int) -> np.ndarray:

        '''
        Parameters
        ----------
            - time : the time from wich read the trajectory

        Returns
        ----------
            - trj : a numpy array of shape (N,2) with the x,y components of the trajectory 
        
        '''

        dump = f"{self.path}/xyz.dump.{str(time)}.{self.extension}"

        # try to read the files
        
        trj =  np.loadtxt(dump, skiprows=9, ndmin=1)       
        trj = trj[:,self.POS_INDICES]
        
        return trj