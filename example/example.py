

import awe

import numpy as np
import os

cfg           = awe.workqueue.Config()
cfg.name      = 'awe-example'
cfg.fastabort = 2
cfg.restarts  = 3

### Don't use this unless MSMBuilder is installed and accessible from the workers.
### Double-check the workerfiles/env.sh script before using execute-task.sh
# cfg.execute ('workerfiles/execute-task.sh' )

### simulate a worker task
cfg.execute ('workerfiles/work.sh' )

### required files
cfg.cache   ('workerfiles/protomol.conf'   )
cfg.cache   ('workerfiles/topol.tpr'       )
cfg.cache   ('workerfiles/with-env'        )
cfg.cache   ('workerfiles/env.sh'          )
cfg.cache   ('workerfiles/Gens.lh5'        )
cfg.cache   ('workerfiles/AtomIndices.dat' )
cfg.cache   ('workerfiles/state0.pdb'      )


iterations = 3
nwalkers   = 2
nstates    = 100

system     = awe.System(topology=awe.PDB('workerfiles/state0.pdb'))

print 'Loading cells and walkers'
srcdir = 'conformations'

for i in xrange(nstates):
    cell = awe.Cell(i, color=np.random.randint(0, 2), weight=np.random.random())

    for j in xrange(nwalkers):

        pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
        pdb     = awe.PDB(pdbpath)
        w       = awe.aweclasses.Walker(pdb.coords)
        cell.add_walker(w)

    system.add_cell(cell)


resample = awe.resample.MultiColor(nwalkers)
resample = awe.resample.SaveWeights(resample)
adaptive = awe.AWE( wqconfig   = cfg,
                    system     = system,
                    iterations = iterations,
                    resample   = resample)

adaptive.run()

print 'Run time:', awe.time.time(), 'seconds'
