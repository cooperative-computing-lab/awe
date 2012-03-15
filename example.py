

import awe

import numpy as np
import os

cfg = awe.workqueue.Config()
cfg.name = 'awe-badi'
cfg.fastabort = 3
cfg.restarts = float('inf')


# cfg.execute('test.exe')

cfg.execute('testinput/execute-task.sh')

cfg.cache('testinput/protomol.conf')
cfg.cache('testinput/topol.tpr')
cfg.cache('testinput/with-env')
cfg.cache('testinput/env.sh')
cfg.cache('testinput/Gens.lh5')
cfg.cache('testinput/AtomIndices.dat')
cfg.cache('testinput/state0.pdb')


iterations = 30
nwalkers   = 10
nstates    = 100

system = awe.System(topology = awe.PDB('testinput/state0.pdb'))
partition = awe.SinkStates()

for color, grp in enumerate( [ [0,39,45,58,66,24,71], [6] ] ):
    partition.add(color, *grp)

print 'Loading cells and walkers'
srcdir = '/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs'
for i in xrange(nstates):


    cell  = awe.Cell(i)
    color = partition.color(cell)
    system.add_cell(cell)


    for j in xrange(nwalkers):

        pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
        pdb     = awe.PDB(pdbpath)
        w       = awe.Walker(start=pdb.coords, assignment=i, color=color, weight=np.random.random())
        system.add_walker(w)



resample = awe.resample.MultiColor(nwalkers, partition)
resample = awe.resample.SaveWeights(resample)
adaptive = awe.AWE( wqconfig   = cfg,
                    system     = system,
                    iterations = iterations,
                    resample   = resample)

adaptive.run()

print 'Run time:', awe.time.time(), 's'
