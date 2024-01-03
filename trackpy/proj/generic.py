import os
from ._projects import Subproject


class Generic(Subproject):
    def __init__(self, genpath):
        self.genpath = genpath

    def filepath(self, sub=None, time=None):
        return self.genpath + super().folder_structure(sub, time)
