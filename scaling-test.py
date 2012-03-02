

import awe
import mdtools

import numpy as np
import os

mdtools.prody.setVerbosity('error')

N_SIM_STEPS = 1000

cfg = awe.workqueue.Config()
cfg.name = 'awe-scaling-test'
cfg.fastabort = 3
cfg.restarts = 0


prefix   = 'scaling-test/dhfr-solv'
inputdir = os.path.join(prefix, 'input')

cfg.execute(os.path.join(inputdir, 'run.sh'))
cfg.cache(os.path.join(inputdir  , '%ssteps.tpr' % N_SIM_STEPS))

iterations = 10
nwalkers   = 10
nstates    = 10

pdb        = mdtools.prody.parsePDB(os.path.join(inputdir, 'conf.pdb'))
walkers    = awe.WalkerGroup(count    = nwalkers * nstates, topology = pdb)

for _ in xrange(nstates * nwalkers):
    w          = awe.Walker( start = pdb.getCoords(),
                             weight = np.random.random(),
                             cell = 0
                             )
    walkers.add(w)


resample = awe.resample.Identity()
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               walkers    = walkers,
                               iterations = iterations,
                               resample   = resample)

adaptive.statsdir = os.path.join(prefix, '%ssteps-long-stats' % N_SIM_STEPS)
adaptive.run()
