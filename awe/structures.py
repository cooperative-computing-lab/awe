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

    """
    awe.structures.PDB

    Reads, manages, and writes PDB files using the ProDy API.

    Fields:
        _pdb - a ProDy AtomGroup object representing the contents of a PDB file

    Methods:
        _from_str    - loads a ProDy AtomGroup from a string representation
        _from_file   - loads a ProDy AtomGroup from a PDB file
        __getstate__ - transform a ProDy AtomGroup into a pickleable state
        __setstate__ - interpret a ProDy AtomGroup from a pickled state
        copy         - make a copy of the ProDy AtomGroup
    """

    def __init__(self, string=None, pdb=None):

        """
        awe.structures.PDB.__init__

        Initialize a new PDB object from a source depending on the parameters.

        Parameters:
            string - a filepath or a string representation of a ProDy AtomGroup
            pdb    - a ProDy AtomGroup

        Returns:
            None
        """

        if string and os.path.exists(string):
            # Load from file if string is a filepath
            self._pdb = PDB._from_file(string)
        elif string:
            # Load from string if string is a string representation
            self._pdb = PDB._from_str(string)
        elif pdb is not None:
            # pdb must be an AtomGroup
            self._pdb = pdb
        else:
            # Otherwise, there is nothing to load
            self._pdb = None

    @staticmethod
    @returns(prody.AtomGroup)
    def _from_str(string):

        """
        awe.structures.PDB._from_str

        Load a ProDy AtomGroup form a string representation.

        Parameters:
            string - a string containing the contents of a ProDy AtomGroup

        Returns:
            A ProDy AtomGroup
        """

        ss = io.StringStream(string)
        return prody.parsePDBStream(ss)

    @staticmethod
    @returns(prody.AtomGroup)
    def _from_file(path):
        
        """
        awe.structures.PDB._from_file

        Load a ProDy AtomGroup from a pdb file.

        Parameters:
            path - a string containing a filepath to a PDB file

        Returns:
            A ProDy AtomGroup
        """

        return prody.parsePDB(path)

    def __getstate__(self):
        """
        See Python docs on the pickle module for more info on __getstate__
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
        See Python docs on the pickle module for more info on __setstate__
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
        
        """
        awe.structures.PDB.copy

        Make a copy of the ProDy AtomGroup.

        Parameters:
            None

        Returns:
            A new instance of PDB
        """
        return PDB(pdb=self._pdb.copy())

