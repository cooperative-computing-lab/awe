

import awe

import numpy as np
import os

cfg = awe.workqueue.Config()
cfg.name = 'awe-badi'
cfg.fastabort = 3
cfg.restarts = 0


cfg.execute('test.exe')

# cfg.execute('testinput/execute-task.sh')
# cfg.cache('testinput/protomol.conf')
# cfg.cache('testinput/topol.tpr')
# cfg.cache('testinput/with-env')
# cfg.cache('testinput/env.sh')
# cfg.cache('testinput/Gens.lh5')
# cfg.cache('testinput/AtomIndices.dat')
# cfg.cache('testinput/state0.pdb')


iterations = 3
nwalkers   = 2
nstates    = 100

system = awe.System(topology = awe.PDB('testinput/state0.pdb'))

print 'Loading cells and walkers'
srcdir = '/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs'
for i in xrange(nstates):

    if i < nstates/2: core = 0
    else            : core = 1

    cell = awe.Cell(i, core=core)
    system.add_cell(cell)

    for j in xrange(nwalkers):

        pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
        pdb     = awe.PDB(pdbpath)
        w       = awe.aweclasses.Walker(start=pdb.coords, assignment=i, color=core, weight=np.random.random())
        system.add_walker(w)



resample = awe.resample.MultiColor(nwalkers)
# resample = awe.resample.SaveWeights(resample)
adaptive = awe.AWE( wqconfig   = cfg,
                    system     = system,
                    iterations = iterations,
                    resample   = resample)

adaptive.run()

print 'Run time:', awe.time.time(), 's'
