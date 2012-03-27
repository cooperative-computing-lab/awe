

import awe

import numpy as np
import os

cfg = awe.workqueue.Config()
cfg.name = 'awe-badi'
cfg.fastabort = 3
cfg.restarts = 0 # float('inf')


# cfg.execute('test.exe')

cfg.execute('testinput/execute-task.sh')

cfg.cache('testinput/sim.mdp')
cfg.cache('testinput/env.sh')
cfg.cache('testinput/Gens.lh5')
cfg.cache('testinput/AtomIndices.dat')
cfg.cache('testinput/state0.pdb')


iterations = 5
nwalkers   = 2
nstates    = 100

system = awe.System(topology = awe.PDB('testinput/state0.pdb'))

partition = awe.SinkStates()
partition.add(0, *range(0,50))
partition.add(1, *range(50,100))


print 'Loading cells and walkers'
srcdir = '/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs'
for i in xrange(nstates):

    if i < nstates / 3:
        cell = awe.Cell(i, core=0)
    elif i > 2 * nstates / 3:
        cell = awe.Cell(i, core=1)
    else:
        cell = awe.Cell(i)

    color = partition.color(cell)
    system.add_cell(cell)


    for j in xrange(nwalkers):

        pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
        pdb     = awe.PDB(pdbpath)
        w       = awe.Walker(start=pdb.coords, assignment=i, color=color, weight=np.random.random())
        system.add_walker(w)



multicolor = awe.resample.MultiColor(nwalkers, partition)
resample = awe.resample.SaveWeights(multicolor)
adaptive = awe.AWE( wqconfig   = cfg,
                    system     = system,
                    iterations = iterations,
                    resample   = resample)

adaptive.run()

multicolor.save_transitions('transitions.dat')

print 'Run time:', awe.time.time(), 's'
