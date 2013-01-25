# -*- mode: Python; indent-tabs-mode: nil -*-  #
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
from collections import defaultdict


_WALKER_ID = 0
_DEFAULT_COLOR = -1
DEFAULT_CORE = -1

class Walker(object):

    """
    Capture the state of a single walker.
    This include the starting and ending coordinates, and it's assignment

    Relevant fields are:

      *start*      : starting coordinates
      *end*        : ending coordinates
      *assignment* : int
    """

    def __init__(self, start=None, end=None, assignment=None, color=_DEFAULT_COLOR, weight=None, wid=None, cellid=None, initid=None):

        assert not (start is None and end is None), 'start = %s, end = %s' % (start, end)

        self._start      = start
        self._end        = end
        self._assignment = assignment
        self._color      = color
        self._weight     = weight
        self._cellid     = cellid

        if wid is None:
            global _WALKER_ID
            self._id     = _WALKER_ID
            _WALKER_ID  += 1
        else:
            self._id     = wid

        self._initid    = initid or self._id

    def __eq__(self, other):
        if not type(self) is type(other):
            return False

        return \
            (self._start     == other._start).all() and \
            self._assignment == other._assignment   and \
            self._color      == other._color        and \
            self._weight     == other._weight


    def restart(self, weight=None, cellid=None):
        assert self._start is not None
        assert self._end   is not None
        assert weight      is not None

        global _WALKER_ID
        wid =  _WALKER_ID
        _WALKER_ID += 1

        cid = cellid or self._cellid

        return Walker(start      = self._end,
                      end        = None,
                      assignment = self._assignment,
                      color      = self._color,
                      weight     = weight,
                      wid        = wid,
                      cellid     = cid,
                      initid     = self._initid)


    @property
    def id(self):         return self._id

    @property
    def cellid(self):     return self._cellid

    @property
    def initid(self):    return self._initid

    @property
    def start(self):      return self._start

    @property
    def end(self):        return self._end

    @end.setter
    def end(self, crds):  self._end = crds

    @property
    def assignment(self): return self._assignment

    @assignment.setter
    def assignment(self, asn):   self._assignment = asn

    @property
    def color(self):      return self._color

    @color.setter
    def color(self, c):   self._color = c

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
    def __init__(self, wqconfig=None, system=None, iterations=-1, resample=None,
                 statsdir = 'stats',
                 checkpointfile='checkpoint', checkpointfreq=1):

        self.statslogger = stats.StatsLogger('stats.log.bz2')
        self.transitionslogger = stats.StatsLogger('cell-transitions.log.bz2')

        self.wq         = workqueue.WorkQueue(wqconfig, statslogger=self.statslogger)
        self.system     = system
        self.iterations = iterations
        self.resample   = resample

        self.iteration  = 0

        self.stats      = stats.AWEStats(logger=self.statslogger)
        self.statsdir   = statsdir

        self.checkpointfile = checkpointfile
        self.checkpointfreq = checkpointfreq


    def checkpoint(self, path):

        tmpfile = path + '.tmp'
        with open(tmpfile, 'wb') as fd:
            pickle.dump(self, fd, pickle.HIGHEST_PROTOCOL)
        shutil.move(tmpfile, path)


    def save_stats(self, dirname):
        if not os.path.exists(dirname):
            print 'Creating directory', dirname
            os.makedirs(dirname)

        awestats = os.path.join(dirname, 'awestats.npy')
        self.stats.save(awestats)
        self.wq.save_stats(dirname)


    def _submit(self):

        for walker in self.system.walkers:
            task = self._new_task(walker)
            self.wq.submit(task)

    @typecheck(Walker)
    @returns(workqueue.WQ.Task)
    def _new_task(self, walker):
        task = self.wq.new_task()
        tag  = self.encode_task_tag(walker)
        task.specify_tag(tag)
        self.marshal_to_task(walker, task)
        return task

    def _try_duplicate_tasks(self):
        i = 0
        while self.wq.can_duplicate_tasks():
            i += 1
            if i > 20: break
            tag    = self.wq.select_tag()
            print time.asctime(), 'Trying to duplicate tag:', tag
            if tag is None: break
            print time.asctime(), 'Duplicating tag', tag
            wid    = self.decode_from_task_tag(tag)['walkerid']
            walker = self.system.walker(wid)
            task   = self._new_task(walker)
            self.wq.submit(task)

    def _recv(self):

        print time.asctime(), 'Recieving tasks'
        system = self.system
        self.stats.time_barrier('start')
        while not self.wq.empty:
            walker = self.wq.recv(self.marshal_from_task)
            system.set_walker(walker)
            self._try_duplicate_tasks()
        self.stats.time_barrier('stop')
        self.wq.clear_tags()
        print system


    def _resample(self):

        self.stats.time_resample('start')
        self.system = self.resample(self.system)
        self.stats.time_resample('stop')
            

    def run(self):
        """
        Run the algorithm
        """

        assert len(self.system.cells  ) > 0
        assert len(self.system.walkers) > 0

        t = time.time()
        self.statslogger.update(t, 'AWE', 'start_unix_time', t)

        try:
            while True:

                if self.iteration >= self.iterations: break

                self.iteration += 1

                print time.asctime(), 'Iteration', self.iteration, 'with', len(self.system.walkers), 'walkers'
                runtime = stats.time.time()
                self.statslogger.update(runtime, 'AWE', 'iteration', self.iteration)
                self.statslogger.update(runtime, 'AWE', 'walkers', len(self.system.walkers))

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

        # finally:
        #     self.save_stats(self.statsdir)


    # @typecheck(int, int)
    # @returns(str)
    def encode_task_tag(self, walker):
        tag = '%(outfile)s|%(cellid)d|%(weight)f|%(walkerid)d' % {
            'outfile' : os.path.join(self.wq.tmpdir, workqueue.RESULT_NAME),
            'cellid'  : walker.assignment,
            'weight'  : walker.weight,
            'walkerid' : walker.id}

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

        

    @typecheck(int, workqueue.WQ.Task)
    def marshal_to_task(self, walker, task):
        

        ### create the pdb
        top        = self.system.topology
        top.coords = walker.start
        pdbdat     = str(top)

        ### send walker to worker
        wdat = pickle.dumps(walker)
        task.specify_buffer(pdbdat, workqueue.WORKER_POSITIONS_NAME, cache=False)
        task.specify_buffer(wdat  , workqueue.WORKER_WALKER_NAME   , cache=False)

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

            walkerstr         = tar.extractfile(workqueue.WORKER_WALKER_NAME).read()
            pdbstring         = tar.extractfile(workqueue.RESULT_POSITIONS).read()
            cellstring        = tar.extractfile(workqueue.RESULT_CELL     ).read()

            pdb               = structures.PDB(pdbstring)
            coords            = pdb.coords
            cellid            = int(cellstring)

            walker            = pickle.loads(walkerstr)

            selftransition = walker.assignment == cellid
            print time.asctime(), 'Walker', walker.id, 'cell', walker.assignment, '->', cellid, selftransition
            self.transitionslogger.update(time.time(), 'AWE', 'cell_transition',
                                          'iteration %s from %s to %s %s' % \
                                              (self.iteration, walker.assignment, cellid, selftransition))

            walker.end        = coords
            walker.assignment = cellid

        os.unlink(result.tag)
        return walker





class Cell(object):

    def __init__(self, cid, weight=1., core=DEFAULT_CORE, walkers=None):
        self._id      = cid
        self._core    = core


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
        return '<System: topology=%s, ncells=%s, nwalkers=%s>' % \
            (type(self.topology), len(self._cells), len(self._walkers))

    def __repr__(self):
        return 'System(topology=%r, cells=%r)' % (self._topology, self._cells)

    def __iadd__(self, other):
        self._cells.update(other._cells)
        self._walkers.update(other._walkers)

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

    @returns(Walker)
    def walker(self, wid):
        return self._walkers[wid]

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

    def clone(self, cells=True):
        _cells = self._cells if cells else dict()
        return System(topology=self.topology, cells=_cells)


class SinkStates(object):

    def __init__(self):
        self._color_state = defaultdict(set)
        self._state_color = dict()

    def add(self, color, *states):
        for state in states:
            self._color_state[color].add(state)
            self._state_color[state] = color

    def color(self, cell):
        if cell.id in self._state_color:
            return self._state_color[cell.id]
        else:
            global _DEFAULT_COLOR
            return _DEFAULT_COLOR

    def states(self, color):
        return self._color_state[color]

    @property
    def ncolors(self): return len(self._color_state)
