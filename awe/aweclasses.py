
import io, stats, workqueue
from util import typecheck


import mdtools

import numpy as np

import os, time


class Walker(object):

    """
    Capture the state of a single walker.
    This include the starting and ending coordinates, the cell, the
    weight, and optional color.

    Relevant fields are:

      *start*  : numpy.ndarray (natoms, 3) : starting coordinates
      *end*    : numpy.ndarray (natoms, 3) : ending coordinates
      *weight* : float
      *color*  : int
      *cell*   : int
    """

    def __init__(self, start=None, end=None, weight=1., color=0, cell=None, wid=-1):

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

    """
    Efficiently store the state between iterations of the AWE
    algorithms.  This class implements the __get/setitem__ methods
    which accept and return Walker instances.

    This is the data structure that is passed to the resampling
    algorithm, which needs to return a new WalkerGroup instance.


    Relevant fields are:

      *startcoords* : numpy.ndarray (count, natoms, 3)
      *endcoords*   : numpy.ndarray (count, natoms, 3)
      *weights*     : numpy.ndarray (count)
      *colors*      : numpy.ndarray (count)
      *cells*       : numpy.ndarray (count)

    """

    @typecheck(count=int, topology=mdtools.prody.AtomGroup)
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


    @typecheck(mdtools.prody.AtomGroup)
    def topology(self, pdb):
        """
        Set the topology
        """

        self.topology = pdb

    @typecheck(Walker)
    def add(self, walker, ix=None):
        """
        Add a Walker to this group
        """

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
        """
        Set a Walker *walker* at the *i*th position
        """

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
        """
        Get the *i*th walker. Returns an instance of Walker.
        """

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
        """
        Use fancy slicing to select a subset of values from this walker group
        Returns a new WalkerGroup
        """

        g             = WalkerGroup(count=len(slice), topology=self.topology)
        g.startcoords = self.startcoords [slice].copy()
        g.endcoords   = self.endcoords   [slice].copy()
        g.weights     = self.weights     [slice].copy()
        g.colors      = self.colors      [slice].copy()
        g.cells       = self.cells       [slice].copy()
        g._count      = len(g.weights)
        g._ix         = g._count

        return g

    def get_pdb(self, k):
        """
        Get a pdb (coordinates + topology) for the *k*th walker
        Returns a mdtools.prody.AtomGroup instance
        """
        pdb = self.topology.copy()
        pdb.setCoords(self.startcoords[k])
        return pdb

    def get_task_params(self, k):
        """
        Get a dictionary of the parameters suitable for creating a WorkQueue Task
        """

        ss  = io.StringStream()
        mdtools.prody.writePDBStream(ss, self.get_pdb(k))

        return {'pdb'    : ss.read()              ,
                'weight' : str( self.weights [k]) ,
                'color'  : str( self.colors  [k]) ,
                'cell'   : str( self.cells   [k]) ,
                'id'     : str(k)                 }




class AWE(object):

    """
    The main driver for the Adaptive Weighted Ensemble algorithm.
    This class manages the marshaling of workers to/from workers,
    updating the current WalkerGroup, and calling the resampleing
    algorithm.

    When constructing an AWE instance, required parameters include:

      *wqconfig*   : an instance of awe.workqueue.Config
      *walkers*    : an instance of awe.aweclasses.WalkerGroup
      *iterations* : number of iterations to run
      *resample*   : the implementation of the resampling algorithm, a subclass of awe.resample.IResampler


    example:
    >>> cfg        = awe.workqueue.Config()
    >>> walkers    = awe.aweclasses.WalkerGroup(...)
    >>> iterations = 42
    >>> resample   = awe.resample.OneColor()
    >>> adaptive   = awe.aweclasses.AWE(wqconfig=cfg, walkers=walkers, iterations=42, resample=resample)
    >>> adaptive.run()
    """

    @typecheck(wqconfig=workqueue.Config, walkers=WalkerGroup, iterations=int)
    def __init__(self, wqconfig=None, walkers=None, iterations=-1, resample=None):

        assert type(wqconfig) is workqueue.Config
        assert type(walkers) is WalkerGroup
        assert type(iterations) is int
        # TODO: assert type(resample) is

        self.wq         = workqueue.WorkQueue(wqconfig)
        self.walkers    = walkers
        self.iterations = iterations
        self.resample   = resample

        self.stats      = stats.AWEStats()
        self.statsdir   = 'stats'



    def save_stats(self, dirname):
        if not os.path.exists(dirname):
            print 'Creating directory', dirname
            os.makedirs(dirname)

        awestats = os.path.join(dirname, 'awestats.npy')
        self.stats.save(awestats)
        self.wq.save_stats(dirname)


    def _submit(self):

        for i in xrange(len(self.walkers)):
            params = self.walkers.get_task_params(i)
            task   = self.wq.new_task(params)
            self.wq.submit(task)


    def _recv(self):

        print time.asctime(), 'Recieving tasks'
        self.stats.time_barrier('start')
        while not self.wq.empty:
            walker                  = self.wq.recv()
            self.walkers[walker.id] = walker
        self.stats.time_barrier('stop')


    def _resample(self):

        self.stats.time_resample('start')
        self.walkers = self.resample(self.walkers)
        self.stats.time_resample('stop')
            

    def run(self):
        """
        Run the algorithm
        """

        try:
            for iteration in xrange(self.iterations):

                print time.asctime(), 'Iteration', iteration

                self.stats.time_iter('start')

                self._submit()
                self._recv()     ## barrier
                self._resample()

                self.stats.time_iter('stop')

        except KeyboardInterrupt:
            pass

        finally:
            self.save_stats(self.statsdir)
