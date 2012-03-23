#!/usr/bin/env python

import tables
import numpy as np
import os, sys

PRECISION = 1000.

hdfile = sys.argv[1]
txtfile = sys.argv[2]


F = tables.openFile(hdfile)
data = F.getNode('/XYZList').read()

### convert from the lossy integer scheme
data = data.astype('float32') / PRECISION

ncells, natoms, dim = data.shape

with open(txtfile, 'w') as fd:
    fd.write(str(ncells) + '\n')
    fd.write(str(natoms) + '\n')
    fd.write(str(dim)    + '\n')

    for c in xrange(ncells):
        for a in xrange(natoms):
            s = ' '.join(map(str, data[c][a]))
            fd.write(s + '\n')

F.close()
