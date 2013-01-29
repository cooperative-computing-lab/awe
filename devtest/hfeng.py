

import awe
import mdtools

import numpy as np
import os

mdtools.prody.setVerbosity('error')


cfg = awe.workqueue.Config()
cfg.fastabort = 2
cfg.execute('testinput/execute-task.sh')
cfg.cache('testinput/protomol.conf')
cfg.cache('testinput/topol.tpr')
cfg.cache('testinput/with-env')
cfg.cache('testinput/env.sh')
cfg.cache('testinput/Gens.lh5')
cfg.cache('testinput/AtomIndices.dat')
cfg.cache('testinput/state0.pdb')


nwalkers = 2
nstates  = 100
walkers  = awe.aweclasses.WalkerGroup(count    = nwalkers * nstates,
                                      topology = mdtools.prody.parsePDB('testinput/state0.pdb'))
weights = [0.062626,0.000196,0.004066,0.003652,0.000015,0.000084,0.014498,0.001696,0.038257,0.003160,0.000106,0.001868,0.001470,0.000137,0.002753,0.011305,0.002285,0.009026,0.000075,0.004512,0.018624,0.000441,0.000088,0.002812,0.040838,0.011943,0.013997,0.004252,0.000101,0.006888,0.004377,0.000812,0.001102,0.006666,0.000098,0.003274,0.022631,0.012045,0.001266,0.037058,0.000249,0.012986,0.011995,0.002893,0.001684,0.059684,0.003408,0.007592,0.002013,0.002243,0.000041,0.004879,0.000009,0.000002,0.018909,0.018078,0.009137,0.000790,0.053345,0.000177,0.004375,0.000184,0.026106,0.024748,0.000249,0.000081,0.059010,0.005670,0.002913,0.003550,0.019244,0.033933,0.002085,0.030817,0.004204,0.006954,0.031398,0.001762,0.008171,0.000423,0.000028,0.001744,0.000037,0.010025,0.037970,0.000008,0.002594,0.019113,0.002565,0.001019,0.000008,0.000064,0.003362,0.026468,0.022342,0.003351,0.005726,0.009862,0.012644,0.013981]

srcdir = '/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs'
for i in xrange(nstates):
    weight = weights[i]
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


resample = awe.resample.OneColor_SaveWeights(nwalkers , datfile = 'weight.dat')
adaptive = awe.aweclasses.AWE( wqconfig   = cfg,
                               walkers    = walkers,
                               iterations = 50,
                               resample   = resample)

adaptive.run()
