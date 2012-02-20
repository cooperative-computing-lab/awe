

import awe
import mdtools

import numpy as np
import os

awe.io.TRACE = False
mdtools.prody.setVerbosity('warning')

cfg = awe.workqueue.Config()
cfg.fastabort = 4
cfg.execute('testinput/execute-task.sh')
cfg.cache('testinput/protomol.conf')
cfg.cache('testinput/topol.tpr')
cfg.cache('testinput/with-env')
cfg.cache('testinput/env.sh')
cfg.cache('testinput/Gens.lh5')
cfg.cache('testinput/AtomIndices.dat')
cfg.cache('testinput/state0.pdb')


nwalkers = 10
nstates  = 100
walkers  = awe.aweclasses.WalkerGroup(count    = nwalkers * nstates,
                                      topology = mdtools.prody.parsePDB('testinput/state0.pdb'))

srcdir = '/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs'
for i in xrange(nstates):
    weight = np.random.random()
    for j in xrange(nwalkers):

        pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
        pdb     = mdtools.prody.parsePDB(pdbpath)
        color   = 0
        cell    = i
        w       = awe.aweclasses.Walker(start  = pdb.getCoords(),
                                        weight = weight,
                                        color  = color,
                                        cell   = cell
                                        )
        walkers.add(w)


resample = awe.resample.OneColor_PlotCellRMSD(nwalkers, ref='testinput/state0.pdb', plotfile='/tmp/test.png')
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               walkers    = walkers,
                               iterations = 10,
                               resample   = resample)

adaptive.run()
resample.plot()
