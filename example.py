# -*- mode: Python; indent-tabs-mode: nil -*-  #

import awe

import numpy as np
import os

cfg = awe.workqueue.Config()
cfg.name = 'test-awe'
cfg.fastabort = 9
cfg.restarts = float('inf')
cfg.maxreps = float('inf')
cfg.debug = 'all'

# cfg.execute('test.exe')

cfg.execute('testinput/execute-task.sh')

cfg.cache('awedata/binaries/$OS-$ARCH/pdb2gmx')
cfg.cache('awedata/binaries/$OS-$ARCH/grompp')
cfg.cache('awedata/binaries/$OS-$ARCH/mdrun')
cfg.cache('awedata/binaries/$OS-$ARCH/assign')

cfg.cache('awedata/gmxtopologies')
cfg.cache('testinput/sim.mdp')
cfg.cache('testinput/env.sh')
cfg.cache('testinput/cells.dat')
cfg.cache('testinput/CellIndices.dat')
cfg.cache('testinput/StructureIndices.dat')
cfg.cache('testinput/state0.pdb')


iterations = 5
nwalkers   = 4
nstates    = 100

weights = np.random.random((nstates,nwalkers))
weights /= np.sum(weights.flatten())

system = awe.System(topology = awe.PDB('testinput/state0.pdb'))

partition = awe.SinkStates()
partition.add(0, *range(0,50))
partition.add(1, *range(50,100))


print 'Loading cells and walkers'
srcdir = 'awedata/pdbs/ala'
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
        w       = awe.Walker(start=pdb.coords, assignment=i, color=color, weight=weights[i,j], cellid=cell.id)
        system.add_walker(w)



multicolor = awe.resample.MultiColor(nwalkers, partition)
resample = awe.resample.SaveWeights(multicolor)
adaptive = awe.AWE( wqconfig   = cfg,
                    system     = system,
                    iterations = iterations,
                    resample   = resample,
                    statsdir   = '/tmp/awe-stats',
                    checkpointfile = '/tmp/awe-checkpoint.dat',
                    checkpointfreq = 1)

adaptive.run()

multicolor.save_transitions('transitions.dat')

print 'Run time:', awe.time.time(), 's'
