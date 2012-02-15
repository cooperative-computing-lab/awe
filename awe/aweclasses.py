
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

        assert type(coords) is np.ndarray, 'Got %s but expected %s' % (type(coords), np.array)
        assert type(weight) is float
        assert type(color)  is Color
        # assert type(cell)   is int
        assert type(wid)    is int

        self.coords = coords
        self.weight = weight
        self.color  = color
        self.cell   = cell
        self.id     = wid

    natoms = property(lambda self: len(self.coords))


class WalkerGroup(object):

    def __init__(self, count=None, topology=None, walkers=list()):

        assert type(count)                                      is int
        assert type(topology)                                   is mdtools.prody.AtomGroup # from a PDB file
        assert len(walkers) == 0 or  type(iter(walkers).next()) is Walker

        self._count    = count

        self.topology  = topology
        self.positions = -1 * np.ones((count, self.natoms, 3))
        self.weights   =      np.ones(count)
        self.colors    =      np.empty(count, dtype=str)
        self.cells     =      np.zeros(count, dtype=int)

        self._ix       = 0

        map(self.add, walkers)


    ### some read-only properties
    natoms = property(lambda self: len(self.topology))


    def topology(self, pdb):

        assert type(pdb) is mdtools.prody.AtomGroup
        self.topology = pdb

    def add(self, walker, ix=None):

        assert type(walker) is Walker

        if ix is None: i = self._ix
        else:          i = ix

        self[i] = walker

        if ix is None:
            self._ix    += 1

        assert self._ix <= self._count

    def __len__(self):
        return self._count

    def __setitem__(self, i, walker):

        assert walker.natoms == self.natoms
        assert i              < self._count

        self.positions[i] = walker.coords
        self.weights  [i] = walker.weight
        self.colors   [i] = str(walker.color)
        self.cells    [i] = walker.cell


    def __getitem__(self, k):
        w = Walker(coords = self.positions    [k]  ,
                   weight = self.weights      [k]  ,
                   color  = Color(self.colors [k]) ,
                   cell   = self.cells        [k]  ,
                   wid    = k                      )
        return w

    def get_pdb(self, k):
        pdb = self.topology.copy()
        pdb.setCoords(self.positions[k])
        return pdb

    def get_task_params(self, k):

        ss  = awe.io.StringStream()
        mdtools.prody.writePDBStream(ss, self.get_pdb(k))
        pdb = ss.read()

        return {'pdb'    : pdb                    ,
                'weight' : str( self.weights [k]) ,
                'color'  :      self.colors  [k]  ,
                'cell'   : str( self.cells   [k]) ,
                'id'     : str(k)                 }




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


    @awe.trace()
    def _submit(self):

        for i in xrange(len(self.walkers)):
            params = self.walkers.get_task_params(i)
            task   = self.wq.new_task(params)
            self.wq.submit(task)


    @awe.trace()
    def _recv(self):

        self.stats.time_barrier('start')
        awe.log('DEBUG: WQ empty?: %s' % self.wq.empty())
        while not self.wq.empty():
            awe.log('Waiting for result')
            walker     = self.wq.recv()
            walkers[k] = walker
        self.stats.time_barrier('stop')


    @awe.trace()
    def _resample(self):

        self.stats.time_resample('start')
        self.walkers = self.resample(self.walkers)
        self.stats.time_resample('stop')
            

    @awe.trace()
    def run(self):

        for iteration in xrange(self.iterations):

            self.stats.time_iter('start')

            self._submit()
            self._recv()     ## barrier
            self._resample()

            self.stats.time_iter('stop')
