"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""

import io, stats, workqueue
from util import typecheck, returns
import structures, util

import numpy as np
import cPickle as pickle

import os, time, shutil


class Walker(object):

    """
    Capture the state of a single walker.
    This include the starting and ending coordinates, and it's assignment

    Relevant fields are:

      *start*      : starting coordinates
      *end*        : ending coordinates
      *assignment* : int
    """

    def __init__(self, start=None, end=None, assignment=None, color=None):

        assert not (start is None and end is None), 'start = %s, end = %s' % (start, end)

        self.start      = start
        self.end        = end
        self.assignment = assignment
        self.color      = color


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

        self.iteration  = 0

        self.stats      = stats.AWEStats()
        self.statsdir   = 'stats'

        self.checkpointfile = 'checkpoint'
        self.checkpointfreq = 1


    def checkpoint(self, path):

        tmpfile = path + '.tmp'
        with open(tmpfile, 'wb') as fd:
            pickle.dump(self, fd)
        shutil.move(tmpfile, path)


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
            while True:

                if self.iteration > self.iterations: break

                self.iteration += 1

                print time.asctime(), 'Iteration', self.iteration

                self.stats.time_iter('start')

                self._submit()
                self._recv()     ## barrier
                self._resample()

                self.stats.time_iter('stop')

                if self.iteration % self.checkpointfreq == 0:
                    print time.asctime(), 'Checkpointing to', self.checkpointfile
                    self.checkpoint(self.checkpointfile)


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

    @typecheck(str)
    @returns(dict)
    def decode_from_task_tag(self, tag):
        split = tag.split('|')
        outfile, cellid, weight, walkerid = tag.split('|')
        return {'cellid'   : int(cellid)   ,
                'weight'   : float(weight) ,
                'walkerid' : int(walkerid) ,
                'outfile'  : outfile       }

        

    @typecheck(int, int, workqueue.WQ.Task)
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

    @typecheck(workqueue.WQ.Task)
    @returns(Walker)
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

    def __init__(self, cid, weight=1., core=0, walkers=None):
        self._id      = cid
        self._weight  = weight
        self._core    = core
        self._walkers = walkers or list()

    @property
    def id(self): return self._id

    @property
    def weight(self): return self._weight

    @property
    def core(self): return self._core

    @property
    def walkers(self): return self._walkers

    @property
    @returns(int)
    def population(self): return len(self._walkers)

    @property
    def color(self, wid):return self._walkers[wid].color


    @typecheck(Walker)
    def add_walker(self, w):
        if w.color is None:
            w.color = self.core
        self._walkers.append(w)

    def walker(self, i):
        return self._walkers[i]

    def as_empty(self, **kws):
        keys = { 'cid'    : kws['cid']    if 'cid'    in kws else self.id    ,
                 'weight' : kws['weight'] if 'weight' in kws else self.weight,
                 'core'   : kws['core']   if 'core'   in kws else self.core  }
        return Cell(**keys)

    def __len__(self):
        return len(self._walkers)

    def __str__(self):
        return '<Cell: %d, weight=%s, core=%s, nwalkers=%s>' % \
            (self.id, self.weight, self.core, len(self.walkers))

    def __repr__(self):
        return str(self)

    def __iter__(self):
        for w in self.walkers:
            yield w



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
    @returns(list)
    def cells(self):
        return self._cells

    @property
    @returns(list)
    def walkers(self):
        ws = list()
        for c in self.cells:
            for w in c:
                ws.append(w)
        return ws

    @property
    # @returns(np.array)
    def weights(self):
        return np.array(map(lambda c: c.weight, self._cells))

    @property
    @returns(set)
    def colors(self):
        colors = set()
        for cell in self.cells:
            for walker in cell:
                colors.add(walker.color)
        return colors


    @typecheck(Cell)
    def add_cell(self, cell):
        if cell.id in set(map(lambda c: c.id, self._cells)):
            raise ValueError, 'Duplicate cell id %d' % cell.id
        self._cells.append(cell)

    @typecheck(Walker)
    def add_walker(self, walker):
        assert walker.assignment >= 0, 'is: %s' % walker.assignment
        self.cell(walker.assignment).add_walker(walker)


    @returns(Cell)
    def cell(self, i):
        cs = filter(lambda c: c.id == i, self._cells)
        assert len(cs) == 1, 'actual: %s' % len(cs)
        return cs[0]

    @returns(bool)
    def has_cell(self, i):
        return len(filter(lambda c: c.id == i, self._cells)) == 1

    # @returns(System)
    def filter_by_cell(self, cellid):
        cells  = filter(lambda c: c.id == cellid, self._cells)
        newsys = System(topology=self._topology, cells=cells)
        return self.clone(cells=cells)

    # @returns(System)
    def filter_by_color(self, color):
        s = self.clone()
        for c in self.cells:
            c2 = c.as_empty()
            for w in c:
                if w.color == color:
                    c2.add_walker(w)
            if len(c2) > 0:
                s.add_cell(c2)
        return s


    # @returns(System)
    def filter_by_core(self, core):
        cells = filter(lambda c: c.core == core, self._cells)
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
