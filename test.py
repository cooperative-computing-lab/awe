

import awe
import mdtools

awe.io.TRACE = False

cfg      = awe.workqueue.Config()
cfg.execute('test.exe')

pdb     = mdtools.prody.parsePDB('topology.pdb')
walkers = awe.aweclasses.WalkerGroup( count    = 1,
                                      topology = pdb )

w       = awe.aweclasses.Walker(      coords = pdb.getCoords(),
                                      weight = 0.,
                                      color  = awe.aweclasses.Color('red'),
                                      cell   = 0 )
walkers.add(w)

resample = awe.resample.DoNothing()
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               walkers    = walkers,
                               iterations = 1,
                               resample   = resample)

adaptive.run()
