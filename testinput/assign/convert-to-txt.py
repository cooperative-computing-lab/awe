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

ncells, ncoords, dim = data.shape

with open(txtfile, 'w') as fd:
    fd.write('ncells: '  + str(ncells) + '\n')
    fd.write('ncoords: ' + str(ncoords) + '\n')
    fd.write('ndims: '   + str(dim)    + '\n')
    fd.write('\n')

    np.savetxt(fd, data.flatten(), fmt='%f')

F.close()
