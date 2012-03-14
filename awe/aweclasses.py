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


_WALKER_ID = 0

class Walker(object):

    """
    Capture the state of a single walker.
    This include the starting and ending coordinates, and it's assignment

    Relevant fields are:

      *start*      : starting coordinates
      *end*        : ending coordinates
      *assignment* : int
    """

    def __init__(self, start=None, end=None, assignment=None, color=None, weight=None, wid=None):

        assert not (start is None and end is None), 'start = %s, end = %s' % (start, end)

        self._start      = start
        self._end        = end
        self._assignment = assignment
        self._color      = color
        self._weight     = weight

        if wid is None:
            global _WALKER_ID
            self._id     = _WALKER_ID
            _WALKER_ID  += 1
        else:
            self._id     = wid


    def restart(self, weight=None):
        assert self._start is not None
        assert self._end   is not None
        assert weight      is not None

        return Walker(start      = self._end,
                      end        = None,
                      assignment = self._assignment,
                      color      = self._color,
                      weight     = self._weight,
                      wid        = self._id)


    @property
    def id(self):         return self._id

    @property
    def start(self):      return self._start

    @property
    def end(self):        return self._end

    @end.setter
    def end(self, coords):
        assert self._end is None
        self._end = coords

    @property
    def assignment(self): return self._assignment

    @assignment.setter
    def assignment(self, asn):   self._assignment = asn

    @property
    def color(self):      return self._color

    @property
    def weight(self):     return self._weight

    @property
    def natoms(self):     return len(self._coords)

    @property
    def ndim(self):       return self._coords.shape[-1]

    @property
    def _coords(self):
        if self.start is not None:
            return self.start
        elif self.end is not None:
            return self.end
        else: raise ValueError, 'Both *start* and *end* should not be None'


    def __str__(self):
        return '<Walker: id=%(id)d, size=%(size)d, dim=%(dim)d, assignment=%(assignment)d, color=%(color)s, weight=%(weight)s>' \
            % {'id' : self.id, 'size'   : self.natoms, 'dim' : self.ndim,
               'assignment' : self.assignment, 'color' : self.color,
               'weight' : self.weight}

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
    def core(self): return self._core

    @property
    def color(self, wid):return self._walkers[wid].color

    def __str__(self):
        return '<Cell: %d, core=%s>' % \
            (self.id, self.core)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if type(other) is not type(self):
            return False

        return \
            self._id     == other._id     and \
            self._core   == other._core



class System(object):

    def __init__(self, topology=None, cells=None):
        self._topology = topology
        self._cells    = cells or dict()
        self._walkers  = dict()


    def __str__(self):
        return '<System: topology=%s, ncells=%s, cells=%s>' % (self.topology, len(self._cells), self._cells)

    def __repr__(self):
        return 'System(topology=%r, cells=%r)' % (self._topology, self._cells)

    def __iadd__(self, other):
        for cell in other.cells:
            if cell.id not in self._cells:
                self.add_cell(cell)
        for walker in other.walkers:
            self.add_walker(walker)

        return self


    @property
    def topology(self): return self._topology.copy()

    @property
    @returns(list)
    def cells(self):
        return self._cells.values()

    @property
    @returns(list)
    def walkers(self):
        return self._walkers.values()

    @property
    # @returns(np.array)
    def weights(self):
        return np.array(map(lambda w: w.weight, self.walkers))

    @property
    @returns(set)
    def colors(self):
        return set(map(lambda w: w.color, self.walkers))

    @typecheck(Cell)
    def add_cell(self, cell):
        if cell.id in self._cells:
            raise ValueError, 'Duplicate cell id %d' % cell.id
        self.set_cell(cell)

    @typecheck(Cell)
    def set_cell(self, cell):
        self._cells[cell.id] = cell

    @typecheck(Walker)
    def add_walker(self, walker):
        assert walker.assignment >= 0, 'is: %s' % walker.assignment
        self.set_walker(walker)

    @typecheck(Walker)
    def set_walker(self, walker):
        self._walkers[walker.id] = walker

    @returns(Cell)
    def cell(self, i):
        return self._cells[i]

    @returns(bool)
    def has_cell(self, i):
        return i in self._cells

    # @returns(System)
    def filter_by_cell(self, cell):
        ws     = filter(lambda w: w.assignment == cell.id, self.walkers)
        newsys = self.clone(cells={cell.id:self.cell(cell.id)})
        for w in ws: newsys.add_walker(w)
        return newsys

    # @returns(System)
    def filter_by_color(self, color):
        ws     = filter(lambda w: w.color == color, self.walkers)
        newsys = self.clone()

        for w in ws:
            newsys.add_walker(w)

            cell = self.cell(w.assignment)
            newsys.set_cell(cell)

        return newsys


    # @returns(System)
    def filter_by_core(self, core):
        cells  = filter(lambda c: c.core == core, self.cells)
        cs     = {}
        for c in cells: cs[c.id] = c

        newsys = self.clone(cells=cs)
        for w in self.walkers:
            if w.assignment in cs:
                newsys.add_walker(w)

        return newsys

    def clone(self, cells=False):
        _cells = self._cells if cells else dict()
        return System(topology=self.topology, cells=_cells)
