
import awe

import mdtools

import numpy as np



class Color(object):
    def __init__(self, name):
        self._name = name
    def __str__(self):
        return self._name
    def __repr__(self):
        return 'Color(%s)' % self._name



class Walker(object):

    def __init__(self, coords=None, weight=0., color=None, cell=None, wid=-1):

        assert type(pdb)    is np.array
        assert type(weight) is float
        assert type(color)  is Color
        assert type(cell)   is int
        assert type(wid)    is int

        self.coords = coords
        self.weight = weight
        self.color  = color
        self.cell   = cell
        self.id     = wid

    natoms = property(lambda self: len(self._coords))


class WalkerGroup(object):

    def __init__(self, count=None, natoms=None, topology=None, walkers=None):

        assert type(count) is int
        assert type(natoms) is int
        assert type(topology) is mdtools.prody.AtomGroup
        # TODO: assert type(walkers) is

        self._count    = count
        self._natoms   = natoms

        self.topology  = topology
        self.positions = -1 * np.ones((count, natoms, 3))
        self.weights   =      np.ones(count)
        self.colors    =      np.empty(count, dtype=str)
        self.cells     =      np.zeros(count, dtype=int)

        self._ix       = 0

        if walkers:
            map(self.add, walkers)


    def topology(self, pdb):

        assert type(pdb) is mdtools.prody.AtomGroup
        self.topology = pdb

    def add(self, walker, ix=None):

        assert type(walker) is Walker

        if ix is None: i = self._ix
        else:          i = ix

        self[i] = walker

        if ix is None:
            self._ix          += 1

    def __setitem__(self, i, walker):

        assert walker.natoms == self._natoms
        assert i              < self._count

        self.positions[ix] = walker.coords
        self.weights  [ix] = walker.weight
        self.colors   [ix] = str(walker.color)
        self.cells    [ix] = walker.cell


    def __getitem__(self, k):
        w = Walker(coords = self.positions    [k],
                   weight = self.weights      [k],
                   color  = Color(self.colors [k]),
                   cell   = self.cells        [k])
        return w

    def get_pdb(self, k):
        pdb = self.topology.copy()
        pdb.setCoords(self.positions[k])
        return pdb

    def get_task_params(self, k):

        ss  = awe.io.StringStream()
        pdb = mdtools.prody.writePDBStream(ss, self.get_pdb(k))

        return {'pdb'    : pdb                    ,
                'weight' : str( self.weights [k]) ,
                'color'  :      self.colors  [k]  ,
                'cell'   : str( self.cells   [k]) }




class AWE(object):

    def __init__(self, wqconfig=None, cells=None, walkers=None, iterations=-1, resample=None):

        assert type(wqconfig) is awe.workqueue.Config
        # TODO: assert type(cells) is
        assert type(walkers) is WalkerGroup
        assert type(iterations) is int
        # TODO: assert type(resample) is

        self.wq         = awe.workqueue.WorkQueue(wqconfig)
        self.cells      = cells
        self.walkers    = walkers
        self.iterations = iterations
        self.resample   = resample

        self.stats      = awe.stats.AWEStats()


    def _submit(self):

        for i in xrange(len(self.walkers)):
            task = self.wq.new_task(self.walkers, i)

            self.wq.submit(task)


    def _recv(self):

        self.stats.time_barrier('start')
        while not self.wq.empty:
            awe.log('Waiting for result')
            walker     = self.wq.recv()
            walkers[k] = walker
        self.stats.time_barrier('stop')


    def _resample(self):

        self.stats.time_resample('start')
        self.walkers = self.resample(self.walkers)
        self.stats.time_resample('stop')
            

    def run(self):

        for iteration in xrange(self.iterations):

            self.stats.time_itr('start')

            self._submit()
            self._recv()     ## barrier
            self._resample()

            self.stats.time_itr('stop')
