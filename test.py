

import awe
import mdtools

import numpy as np

awe.io.TRACE = False

cfg      = awe.workqueue.Config()
cfg.execute('test.exe2')

pdb     = mdtools.prody.parsePDB('topology.pdb')
walkers = awe.aweclasses.WalkerGroup( count    = 2,
                                      topology = pdb )

w       = awe.aweclasses.Walker(      coords = pdb.getCoords(),
                                      weight = np.random.random(),
                                      color  = 0,
                                      cell   = 0 )
walkers.add(w)

w       = awe.aweclasses.Walker(      coords = pdb.getCoords(),
                                      weight = np.random.random(),
                                      color  = 0,
                                      cell   = 1 )
walkers.add(w)

resample = awe.resample.Simple(2)
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               walkers    = walkers,
                               iterations = 2,
                               resample   = resample)

adaptive.run()
