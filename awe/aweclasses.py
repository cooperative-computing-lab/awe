
import io, stats, workqueue
from util import typecheck, returns
import structures

import numpy as np

import os, time


class Walker(object):

    """
    Capture the state of a single walker.
    This include the starting and ending coordinates, and it's assignment

    Relevant fields are:

      *start*      : starting coordinates
      *end*        : ending coordinates
      *assignment* : int
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
        return '<Walker: assignment=%(assignment)d>' % {'assignment' : self.assignment}

    def __repr__(self): return str(self)




class AWE(object):

    """
    The main driver for the Adaptive Weighted Ensemble algorithm.
    This class manages the marshaling of workers to/from workers,
    updating the current WalkerGroup, and calling the resampleing
    algorithm.

    When constructing an AWE instance, required parameters include:

      *wqconfig*   : an instance of awe.workqueue.Config
      *system*    : an instance of awe.aweclasses.System
      *iterations* : number of iterations to run
      *resample*   : the implementation of the resampling algorithm, a subclass of awe.resample.IResampler

    Some relevant fields

      *stats*    : the instance of awe.stats.AWEStats
      *statsdir* : location to save runtime statistics
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

        finally:
            self.save_stats(self.statsdir)


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
        top        = self.system.topology
        top.coords = walker.start
        pdbdat     = str(top)

        ### send walker to worker
        task.specify_buffer(pdbdat, workqueue.WORKER_POSITIONS_NAME, cache=False)

        ### specify output
        self.specify_task_output_file(task)

    def specify_task_output_file(self, task):
        output = os.path.join(self.wq.tmpdir, task.tag)
        task.specify_output_file(output, remote_name = workqueue.WORKER_RESULTS_NAME, cache=False)

    # @typecheck(workqueue.WQ.Task)
    # @returns(Walker)
    def marshal_from_task(self, result):

        import tarfile
        with tarfile.open(result.tag) as tar:

            pdbstring  = tar.extractfile(workqueue.RESULT_POSITIONS).read()
            cellstring = tar.extractfile(workqueue.RESULT_CELL     ).read()

            pdb    = structures.PDB(pdbstring)
            coords = pdb.coords
            cellid = int(cellstring)

            walker = Walker(end=coords, assignment=cellid)

        os.unlink(result.tag)
        return walker





class Cell(object):

    def __init__(self, cid, weight=1., color=0, walkers=None):
        self._id      = cid
        self._weight  = weight
        self._color   = color
        self._walkers = walkers or list()

    @property
    def id(self): return self._id

    @property
    def weight(self): return self._weight

    @property
    def color(self): return self._color

    @property
    def walkers(self): return self._walkers

    @property
    def population(self): return len(self._walkers)

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

    def __repr__(self):
        return str(self)



class System(object):

    def __init__(self, topology=None, cells=None):
        self._topology = topology
        self._cells    = cells or list()


    def __str__(self):
        return '<System: topology=%s, ncells=%s, cells=%s>' % (self.topology, len(self._cells), self._cells)

    def __repr__(self):
        return 'System(topology=%r, cells=%r)' % (self._topology, self._cells)

    def __iadd__(self, other):
        for cell in other.cells:
            if not self.has_cell(cell.id):
                self.add_cell(cell)
            else:
                mycell = self.cell(cell.id)
                for walker in cell.walkers:
                    mycell.add_walker(walker)

        return self


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
        assert walker.assignment >= 0, 'is: %s' % walker.assignment
        self.cell(walker.assignment).add_walker(walker)

    # @typecheck(int)
    # @returns(Cell)
    def cell(self, i):
        cs = filter(lambda c: c.id == i, self._cells)
        assert len(cs) == 1, 'actual: %s' % len(cs)
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
    def clone(self, cells=None):
        cells = cells or list()
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
