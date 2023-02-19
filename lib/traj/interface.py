from abc import ABC, abstractmethod

class ReadInterface(ABC):
    ''' Abstract interface to read clusters from file'''
    @abstractmethod
    def __init__(self, path, extension=None) -> None:
        pass

    @abstractmethod
    def read(self, time):
        ''' Abstract method to read simulation.
        It expects in output a list of particles'''
        pass
    
