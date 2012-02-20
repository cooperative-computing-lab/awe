
import awe

import mdtools

import numpy as np

import time


class Walker(object):

    def __init__(self, start=None, end=None, weight=0., color=None, cell=None, wid=-1):

        assert not (start is None and end is None)

        self.start  = start
        self.end    = end
        self.weight = weight
        self.color  = color
        self.cell   = cell
        self.id     = wid

    @property
    def natoms(self):
        if self.start is not None:
            return len(self.start)
        elif self.end is not None:
            return len(self.end)
        else: raise ValueError, 'Both *start* and *end* should not be None'

    def __str__(self):
        return \
            'Walker : weight=%(weight)r, color=%(color)r, cell=%(cell)r, wid=%(wid)r)' % \
            {'weight' : self.weight, 'color' : self.color, 'cell' : self.cell, 'wid' : self.id}


    def __repr__(self):
        return \
            'Walker(start=%(start)r, end=%(end)r, weight=%(weight)r, color=%(color)r, cell=%(cell)r, wid=%(wid)r)' \
            % {'start' : self.start, 'end' : self.end, 'weight' : self.weight, 'color' : self.color,
               'cell' : self.cell, 'wid' : self.id}


class WalkerGroup(object):

    @awe.typecheck(count=int, topology=mdtools.prody.AtomGroup)
    def __init__(self, count=None, topology=None, walkers=list()):

        assert len(walkers) == 0 or  type(iter(walkers).next()) is Walker

        self._count      = count

        self.topology    = topology
        self.startcoords = np.zeros((count, self.natoms, 3))
        self.endcoords   = np.zeros((count, self.natoms, 3))
        self.weights     = np.ones(count)
        self.colors      = np.zeros(count, dtype=int)
        self.cells       = np.zeros(count, dtype=int)

        self._ix         = 0

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

    def __iter__(self):
        for i in xrange(self._ix):
            yield self[i]

    def __setitem__(self, i, walker):

        assert not (walker.start is None and walker.end is None)
        assert walker.natoms == self.natoms
        assert i              < self._count

        if walker.start is not None:
            self.startcoords[i] = walker.start

        if walker.end is not None:
            self.endcoords[i]   = walker.end

        self.weights  [i]       = walker.weight
        self.colors   [i]       = walker.color
        self.cells    [i]       = walker.cell


    def __getitem__(self, k):
        start = None if (self.startcoords[k] == 0).all() else self.startcoords[k]
        end   = None if (self.endcoords  [k] == 0).all() else self.endcoords  [k]

        w = Walker(start  = start            ,
                   end    = end              ,
                   weight = self.weights [k] ,
                   color  = self.colors  [k] ,
                   cell   = self.cells   [k] ,
                   wid    = k                )
        return w

    def getslice(self, slice):

        g             = WalkerGroup(count=len(slice), topology=self.topology)
        g.startcoords = self.startcoords [slice]
        g.endcoords   = self.endcoords   [slice]
        g.weights     = self.weights     [slice]
        g.colors      = self.colors      [slice]
        g.cells       = self.cells       [slice]
        g._count      = len(g.weights)
        g._ix         = g._count

        return g

    def get_pdb(self, k):
        pdb = self.topology.copy()
        pdb.setCoords(self.startcoords[k])
        return pdb

    def get_task_params(self, k):

        ss  = awe.io.StringStream()
        mdtools.prody.writePDBStream(ss, self.get_pdb(k))

        return {'pdb'    : ss.read()              ,
                'weight' : str( self.weights [k]) ,
                'color'  : str( self.colors  [k]) ,
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

        print time.asctime(), 'Recieving tasks'
        self.stats.time_barrier('start')
        while not self.wq.empty:
            walker                  = self.wq.recv()
            self.walkers[walker.id] = walker
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
