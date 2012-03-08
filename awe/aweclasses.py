
import io, stats, workqueue
from util import typecheck, returns


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

    def __init__(self, start=None, end=None, assignment=0):

        assert not (start is None and end is None)

        self.start      = start
        self.end        = end
        self.assignment = assignment

    @property
    def natoms(self):
        return len(self._coords)

    @property
    def ndim(self):
        return self._coords.shape[-1]

    @property
    def _coords(self):
        if self.start is not None:
            return self.start
        elif self.end is not None:
            return self.end
        else: raise ValueError, 'Both *start* and *end* should not be None'


    def __str__(self):
        return 'Walker : assignment=%(assignment)d' % {'assignment' : self.assignment}

    def __repr__(self): return str(self)


    # def __repr__(self):
    #     return \
    #         'Walker(start=%(start)r, end=%(end)r, assignment=%(assignment)r)' \
    #         % {'start' : self.start, 'end' : self.end, 'assignment' : self.assignment}



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

    # @typecheck(count=int, topology=mdtools.prody.AtomGroup)
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


    # @typecheck(mdtools.prody.AtomGroup)
    def topology(self, pdb):
        """
        Set the topology
        """

        self.topology = pdb

    # @typecheck(Walker)
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

    # @typecheck(wqconfig=workqueue.Config, system=System, iterations=int)
    def __init__(self, wqconfig=None, system=None, iterations=-1, resample=None):

        self.wq         = workqueue.WorkQueue(wqconfig)
        self.system     = system
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

        for cell in self.system.cells:
            for wid, walker in enumerate(cell.walkers):
                task = self.wq.new_task()
                task.specify_tag(self.encode_task_tag(cell.id, wid))
                self.marshal_to_task(cell.id, wid, task)
                self.wq.submit(task)



    def _recv(self):

        print time.asctime(), 'Recieving tasks'
        system = self.system.as_empty()
        self.stats.time_barrier('start')
        while not self.wq.empty:
            walker = self.wq.recv(self.marshal_from_task)
            system.add_walker(walker)
        self.stats.time_barrier('stop')


    def _resample(self):

        self.stats.time_resample('start')
        self.walkers = self.resample(self.system)
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

        # finally:
        #     self.save_stats(self.statsdir)


    # @typecheck(int, int)
    # @returns(str)
    def encode_task_tag(self, cid, wid):
        cell = self.system.cell(cid)
        tag = '%(outfile)s|%(cellid)d|%(weight)f|%(walkerid)d' % {
            'outfile' : os.path.join(self.wq.tmpdir, workqueue.RESULT_NAME),
            'cellid'  : cid,
            'weight'  : cell.weight,
            'walkerid' : wid}

        return tag

    # @typecheck(str)
    # @returns(dict)
    def decode_from_task_tag(self, tag):
        split = tag.split('|')
        outfile, cellid, weight, walkerid = tag.split('|')
        return {'cellid'   : int(cellid)   ,
                'weight'   : float(weight) ,
                'walkerid' : int(walkerid) ,
                'outfile'  : outfile       }

        

    # @typecheck(int, int, workqueue.WQ.Task)
    def marshal_to_task(self, cid, wid, task):
        
        cell   = self.system.cell(cid)
        walker = cell.walker(wid)

        ### create the pdb
        ss     = io.StringStream()
        top    = self.system.topology
        if walker.start is None:
            print 'walker.start is None', cid, wid # , walker.end
        top.setCoords(walker.start)
        mdtools.prody.writePDBStream(ss, top)
        pdbdat = ss.read()

        ### send walker to worker
        task.specify_buffer(pdbdat, workqueue.WORKER_POSITIONS_NAME, cache=False)

        ### specify output
        self.specify_task_output_file(task)

    def specify_task_output_file(self, task):
        output = os.path.join(self.wq.tmpdir, workqueue.RESULT_NAME)
        task.specify_output_file(output, remote_name = workqueue.WORKER_RESULTS_NAME, cache=False)

    # @typecheck(workqueue.WQ.Task)
    # @returns(Walker)
    def marshal_from_task(self, result):
        d = self.decode_from_task_tag(result.tag)

        import tarfile
        with tarfile.open(d['outfile']) as tar:

            pdbstring  = tar.extractfile(workqueue.RESULT_POSITIONS).read()
            cellstring = tar.extractfile(workqueue.RESULT_CELL     ).read()

            ss     = io.StringStream(pdbstring)
            pdb    = mdtools.prody.parsePDBStream(ss)
            coords = pdb.getCoords()
            cellid = int(cellstring)

            walker = Walker(end=coords, assignment=cellid)

        # os.unlink(d['outfile'])
        return walker





class Cell(object):

    def __init__(self, cid, weight=1., color=0, walkers=list()):
        self._id      = cid
        self._weight  = weight
        self._color   = color
        self._walkers = walkers

    @property
    def id(self): return self._id

    @property
    def weight(self): return self._weight

    @property
    def color(self): return self._color

    @property
    def walkers(self): return self._walkers

    # @typecheck(Walker)
    def add_walker(self, w):
        self._walkers.append(w)

    def walker(self, i):
        return self._walkers[i]

    def as_empty(self):
        return Cell(self._id, weight=self._weight, color=self._color)

    def __len__(self):
        return len(self._walkers)

    def __str__(self):
        return '<Cell: %d, weight=%s, color=%s, nwalkers=%s>' % \
            (self.id, self.weight, self.color, len(self.walkers))



class System(object):

    def __init__(self, topology=None, cells=list()):
        self._topology = topology
        self._cells    = cells


    def __str__(self):
        return '<System: topology=%s, ncells=%s, cells=%s>' % (self.topology, len(self._cells), self._cells[1])

    def __repr__(self):
        return 'System(topology=%r, cells=%r)' % (self._topology, self._cells)


    @property
    def topology(self): return self._topology.copy()

    @property
    # @returns(list)
    def cells(self):
        return self._cells

    @property
    # @returns(np.array)
    def weights(self):
        return np.array(map(lambda c: c.weight, self._cells))

    @property
    # @returns(set)
    def colors(self):
        return set(map(lambda c: c.color, self._cells))



    # @typecheck(Cell)
    def add_cell(self, cell):
        if cell.id in set(map(lambda c: c.id, self._cells)):
            raise ValueError, 'Duplicate cell id %d' % cell.id
        self._cells.append(cell)

    # @typecheck(Walker)
    def add_walker(self, walker):
        self.cell(walker.assignment).add_walker(walker)

    # @typecheck(int)
    # @returns(Cell)
    def cell(self, i):
        cs = filter(lambda c: c.id == i, self._cells)
        assert len(cs) == 1
        return cs[0]

    # @typecheck(int)
    # @returns(bool)
    def has_cell(self, i):
        return len(filter(lambda c: c.id == i, self._cells)) == 1

    # @typecheck(int)
    # @returns(System)
    def filter_by_cell(self, cellid):
        cells  = filter(lambda c: c.id == cellid, self._cells)
        newsys = System(topology=self._topology, cells=cells)
        return self.clone(cells=cells)

    # @typecheck(int)
    # @returns(System)
    def filter_by_color(self, color):
        cells = filter(lambda c: c.color == color, self._cells)
        return self.clone(cells=cells)

    # @returns(list)
    def empty_cells(self):
        cells = list()
        for cell in self.cells:
            cells.append(cell.as_empty())
        return cells

    # @returns(System)
    def clone(self, cells=list()):
        return System(topology=self.topology, cells=cells)

    # @returns(System)
    def as_empty(self):
        return self.clone(cells=self.empty_cells())



    def as_state(self):

        ### sanity check
        nwalkers = len(self._cells[0])
        natoms   = self._cells[0].walkers[0].natoms
        ndim     = self._cells[0].walkers[0].ndim
        for cell in self._cells:
            assert len(cell) == nwalkers
            for w in cell.walkers:
                assert natoms == w.natoms
                assert ndim   == w.ndim

        ncells      = len(self._cells)
        size        = ncells * nwalkers
        cells       = np.zeros(size, dtype=int)
        weights     = np.ones(size)
        colors      = np.zeros(size, dtype=int)
        startcoords = np.zeros((size, natoms, ndim))
        endcoords   = np.zeros((size, natoms, ndim))
        assignments = np.zeros(size, dtype=int)
        walkers     = np.arange(size)

        ix = 0
        for cell in self._cells:
            for walker in cell.walkers:

                cells       [ix] = cell.id
                weights     [ix] = cell.weight
                colors      [ix] = cell.color
                startcoords [ix] = walker.start
                endcoords   [ix] = walker.end
                assignments [ix] = walker.assignment

                ix += 1

        return State(cells=cells, weights=weights, colors=colors,
                     startcoords=startcoords, endcoords=endcoords, assignments=assignments,
                     walkers=walkers)

class State(object):
    def __init__(self,
                 cells=None, weights=None, colors=None,
                 startcoords=None, endcoords=None, assignments=None,
                 walkers=None):

        self._cells       = cells
        self._weights     = weights
        self._colors      = colors
        self._startcoords = startcoords
        self._endcoords   = endcoords
        self._assignments = assignments
        self._walkers     = walkers

    @property
    def cells(self)       : return self._cells       .copy()

    @property
    def weights(self)     : return self._weights     .copy()

    @property
    def colors(self)      : return self._colors      .copy()

    @property
    def startcoords(self) : return self._startcoords .copy()

    @property
    def endcoords(self)   : return self._endcoords   .copy()

    @property
    def assignments(self) : return self._assignments .copy()

    @property
    def walkers(self)     : return self._walkers     .copy()


    def slice_by(self, test):
        ixs = np.where(test)

        return State(cells       = self.cells[ixs],
                     weights     = self.weights[ixs],
                     colors      = self.colors[ixs],
                     startcoords = self.startcoords[ixs],
                     endcoords   = self.endcoords[ixs],
                     assignments = self.assignments[ixs],
                     walkers     = self.walkers[ixs])

    def walker(self, i):
        return Walker(start=self.startcoords[i], end=self.endcoords[i], assignment=self.assignments[i])
