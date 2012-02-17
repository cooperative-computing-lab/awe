

import awe
import mdtools

import numpy as np

awe.io.TRACE = False

cfg = awe.workqueue.Config()
cfg.execute('testinput/execute-task.sh')
cfg.cache('testinput/protomol.conf')
cfg.cache('testinput/topol.tpr')
cfg.cache('testinput/with-env')
cfg.cache('testinput/env.sh')
cfg.cache('testinput/Gens.lh5')
cfg.cache('testinput/AtomIndices.dat')
cfg.cache('testinput/state0.pdb')

nwalkers = 910
walkers  = awe.aweclasses.WalkerGroup( count    = nwalkers,
                                      topology = mdtools.prody.parsePDB('testinput/state0.pdb'))

for i in xrange(nwalkers):
    pdb    = mdtools.prody.parsePDB('teststartpdbs/state%d.pdb'  % i)
    weight = float(open(      'teststartpdbs/weight%d.txt' % i).read().strip())
    color  = 0
    cell   = i
    w      = awe.aweclasses.Walker(coords = pdb.getCoords(),
                                   weight = weight,
                                   color  = color,
                                   cell   = cell
                                   )
    walkers.add(w)


resample = awe.resample.Simple(20)
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               walkers    = walkers,
                               iterations = 2,
                               resample   = resample)

adaptive.run()
