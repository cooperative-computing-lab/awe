
from util import typecheck, returns
import io

import prody

import numpy as np

import os

class PDB(object):

    def __init__(self, string=None, pdb=None):
        if string and os.path.exists(string):
            self._pdb = prody.parsePDB(string)
        elif string:
            ss = io.StringStream(string)
            self._pdb = prody.parsePDBStream(ss)
        elif pdb is not None:
            self._pdb = pdb
        else:
            self._pdb = None

    def __str__(self):
        ss = io.StringStream()
        prody.writePDBStream(ss, self._pdb)
        return ss.read()

    @property
    @returns(np.ndarray)
    def coords(self):
        return self._pdb.getCoords()

    @coords.setter
    @typecheck(np.ndarray)
    def coords(self, xyz):
        self._pdb.setCoords(xyz)

    def copy(self):
        return PDB(pdb=self._pdb.copy())

