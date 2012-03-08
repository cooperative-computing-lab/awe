

import awe
import mdtools

import numpy as np
import os

mdtools.prody.setVerbosity('error')

timer = awe.stats.Timer()

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


iterations = 30
nwalkers = 2
nstates  = 5
walkers  = awe.aweclasses.WalkerGroup(count    = nwalkers * nstates,
                                      topology = mdtools.prody.parsePDB('testinput/state0.pdb'))

system = awe.System(topology = mdtools.prody.parsePDB('testinput/state0.pdb'))

srcdir = '/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs'
timer.start()
for i in xrange(nstates):
    cell = awe.Cell(i, color=np.random.randint(0, 2), weight=np.random.random())

    for j in xrange(nwalkers):

        pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
        pdb     = mdtools.prody.parsePDB(pdbpath)
        w       = awe.aweclasses.Walker(pdb.getCoords())
        cell.add_walker(w)

    system.add_cell(cell)

timer.stop()

print 'Initialization overhead:', timer.elapsed(), 's'


resample = awe.resample.MultiColor(nwalkers)
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               system     = system,
                               iterations = iterations,
                               resample   = resample)

timer.start()
adaptive.run()
timer.stop()

print 'Run time:', timer.elapsed(), 's'
print 'Total time:', awe.time.time(), 's'
