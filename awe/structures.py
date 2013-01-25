# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


from util import typecheck, returns
import io

import prody

import numpy as np

import os

class PDB(object):

    def __init__(self, string=None, pdb=None):
        if string and os.path.exists(string):
            self._pdb = PDB._from_file(string)
        elif string:
            self._pdb = PDB._from_str(string)
        elif pdb is not None:
            self._pdb = pdb
        else:
            self._pdb = None

    @staticmethod
    @returns(prody.AtomGroup)
    def _from_str(string):
        ss = io.StringStream(string)
        return prody.parsePDBStream(ss)

    @staticmethod
    @returns(prody.AtomGroup)
    def _from_file(path):
        return prody.parsePDB(path)

    def __getstate__(self):
        """
        ProDy's AtomGroup cannot be pickled.
        The workaround is to marshal it to/from a string representation.
        Ugly, but it works.
        """
        pdbstring = str(self)
        odict = self.__dict__.copy()
        odict['_pdb'] = pdbstring
        return odict

    def __setstate__(self, odict):
        """
        Marshal the ProDy AtomGroup in from a string representation.
        See __getstate__ (above) for rationale
        """
        pdbstring = odict['_pdb']
        odict['_pdb'] = PDB._from_str(pdbstring)
        self.__dict__.update(odict)

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

