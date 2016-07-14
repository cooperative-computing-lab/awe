# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""

from . import io_tools, stats, workqueue
from .util import typecheck, returns
from . import structures, util

import trax

import numpy as np
import pickle

import os, time, shutil, random
import tempfile
from collections import defaultdict
import ctypes
import sys
import resource

_WALKER_ID = 0
_DEFAULT_COLOR = -1
DEFAULT_CORE = -1

class Walker(object):

    """
    Container for information about a single walker.

    Fields:
        id         - the id of the walker
        cellid     - the id of the cell that the walker is currently in
        initid     - the id of the walker at initialization
        start      - starting atomic coordinates
        end        - ending atomic coordinates
        assignment - the cell assignment of the walker
        color      - the color of the walker
        natoms     - the number of atoms in the walker
        ndim       - the number of coordinate dimensions

    Methods:
        restart - reset the walker to some initial conditions
    """

    def __init__(self, start=None, end=None, assignment=None, color=_DEFAULT_COLOR, weight=None, wid=None, cellid=None, initid=None):

        """
        Initialize a new instance of Walker.

        Parameters:
            start      - a list of starting atomic coordinates
            end        - a list of ending atomic coordinates
            assignment - the cell to which the walker is assigned
            color      - the color of the walker
            weight     - the weight of the walker
            wid        - the id of the walker
            cellid     - the id if the cell the walker is in
            initid     - the id of the walker at initialization

        Returns:
            None
        """

        # At least one set of atomic coordinates is needed to use the walker
        assert not (start is None and end is None), 'start = %s, end = %s' % \
         (start, end)

        self._start      = start
        self._end        = end
        self._assignment = assignment
        self._color      = color
        self._weight     = weight
        self._cellid     = cellid
        self._valid      = True

        # Assign the walker the next global id number if none was suuplies
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

        """
        Reset a walker to its initial conditions.

        Parameters:
            weight - the new weight to assign the walker
            cellid - the id of the cell the walker is currently in

        Returns:
            A new instance of Walker with the initial conditions of this Walker
        """

        # The walker must have been processed to be reset
        assert self._start is not None
        assert self._end   is not None
        assert weight      is not None

        # Assign the walker a new id
        global _WALKER_ID
        wid =  _WALKER_ID
        _WALKER_ID += 1

        # Tell the walker which cell it is in
        cid = cellid or self._cellid

        # Initialize a new walker with the reset settings
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
    def valid(self):      return self._valid;

    def mark_invalid(self):
        self._valid = False

    @property
    def _coords(self):
        if self.start is not None:
            return self.start
        elif self.end is not None:
            return self.end
        else: raise ValueError('Both *start* and *end* should not be None')


    def __str__(self):
        return '<Walker: id=%(id)d, size=%(size)d, dim=%(dim)d, assignment=%(assignment)d, color=%(color)s, weight=%(weight)s>' \
            % {'id' : self.id, 'size'   : self.natoms, 'dim' : self.ndim,
               'assignment' : self.assignment, 'color' : self.color,
               'weight' : self.weight}

    def __repr__(self): return str(self)




class AWE(object):

    """
    The main driver for the Accelerated Weighted Ensemble algorithm.

    This class manages the marshaling of workers to/from workers,
    updating the current WalkerGroup, and calling the resampleing
    algorithm.

    Fields:
        wq             -
        system         -
        iterations     -
        resample       -
        iteration      -
        currenttask    -
        stats          -
        traxlogger     -
        checkpointfreq -

    Methods:
        checkpoint               -
        logwalker                -
        recover                  -
        run                      -
        encode_task_tag          -
        decode_from_task_tag     -
        marshal_to_task          -
        specify_task_output_file -
        marshal_from_task        -
    """

    # @typecheck(wqconfig=workqueue.Config, system=System, iterations=int)
    def __init__(self, wqconfig=None, system=None, iterations=-1, resample=None,
                 traxlogger = None, checkpointfreq=1, verbose=False, log_it=False):
        """
        Initialize a new instance of AWE.

        Parameters:
            wqconfig       - workqueue.Config instance for WorkQueue settings
            system         - the aweclasses.System instance to run
            iterations     - the number of iterations to run
            resample       - the resampler (typically resample.MultiColor)
            traxlogger     - trax logging utility (see trax module)
            checkpointfreq - how often a restart checkpoint should be recorded

        Returns:
            None
        """
        self._verbose = verbose
        self._print_start_screen()

        self.statslogger = stats.StatsLogger('debug/task_stats.log.gz')
        self.transitionslogger = stats.StatsLogger(
            'debug/cell-transitions.log.gz')

        #self.tmpdir = tempfile.mkdtemp(prefix="awe-tmp.")
        self.wq = workqueue.WorkQueue(wqconfig, statslogger=self.statslogger, log_it=log_it)
        #self.wq.tmpdir = self.tmpdir
        self.system = system
        self.iterations = iterations
        self.resample = resample

        self.iteration = 0

        self.currenttask = 0

        self.stats = stats.AWEStats(logger=self.statslogger)

        self.traxlogger = traxlogger or trax.SimpleTransactional(
                                            checkpoint='debug/trax.cpt',
                                            log='debug/trax.log')

        self.checkpointfreq = checkpointfreq

        self._firstrun  = True
        self._log = log_it

    def _print_start_screen(self):
        """
        Print the software development credits.

        Parameters:
            None

        Returns:
            None
        """
        if self._verbose:
            start_str = "********************************* AWE - Acclerated Weighted Ensemble *******************************\n"
            start_str += "AUTHORS:\n"
            start_str += "  Badi' Abdul-Wahid\n"
            start_str += "  Haoyun Feng\n"
            start_str += "  Jesus Izaguirre\n"
            start_str += "  Eric Darve\n"
            start_str += "  Ronan Costaouec\n"
            start_str += "  Dinesh Rajan\n"
            start_str += "  Douglas Thain\n"
            start_str += "\n"
            start_str += "CITATION:\n"
            start_str += "  Badi Abdul-Wahid, Li Yu, Dinesh Rajan, Haoyun Feng, Eric Darve, Douglas Thain, Jesus A. Izaguirre,\n"
            start_str += "  Folding Proteins at 500 ns/hour with Work Queue,\n"
            start_str += "  8th IEEE International Conference on eScience (eScience 2012), October, 2012.\n"
            start_str += "\n"
            start_str += "WEB PAGE:\n"
            start_str += "  www.nd.edu/~ccl/software/awe\n"
            start_str += "***************************************************************************************************\n"

            print(start_str)

    def checkpoint(self):
        """
        Create a checkpoint from which to restart. Replacement for pickling.

        Parameters:
            None

        Returns:
            None
        """

        cpt = self.traxlogger.cpt_path

        if os.path.exists(cpt):
            shutil.move(cpt, cpt + '.last')

        chk = dict(system         = self.system,
                   iterations     = self.iterations,
                   iteration      = self.iteration,
                   resample       = self.resample,
                   checkpointfreq = self.checkpointfreq
                   )

        chk['_firstrun'] = self._firstrun
        self.traxlogger.checkpoint(chk)


    def logwalker(self, walker):
        """
        Output log information for a walker.

        Parameters:
            walker - a walker to log

        Returns:
            None
        """

        self.traxlogger.log(walker)

    def _trax_log_recover(self, obj, value):
        """
        Recover a walker from its log information. See trax module.

        Parameters:
            obj   - the AWE instance
            value - the walker to recover

        Returns:
            None
        """
        if self._verbose:
           print('Recovering walker', value.id)

        # Add the walker back into the system
        obj['system'].set_walker(value)

    def recover(self):
        """
        Recover from interruption using the last checkpoint. Populates the AWE
        instance using values from the checkpoint.

        Parameters:
            None

        Returns:
            None
        """

        cpt = self.traxlogger.cpt_path

        if os.path.exists(cpt):
            if self._verbose:
                print('Recovering', cpt)

            # Get all attributes from the checkpoint
            parms = self.traxlogger.recover(self._trax_log_recover)

            # Reset all AWE parameters
            for a in parms.keys():
                setattr(self, a, parms[a])


    def _submit(self):
        """
        Send each walker in the system to WorkQueue as a task to be executed.

        Parameters:
            None

        Returns:
            None
        """

        for walker in self.system.walkers:
            if walker.end is None:
                task = self._new_task(walker)
                self.wq.submit(task)

    @typecheck(Walker)
    @returns(workqueue.WQ.Task)
    def _new_task(self, walker):
        """
        Create a new WorkQueue Task and give it the correct files to run.

        Parameters:
            walker - the walker to assign to the task

        Returns:
            A WorkQueue Task instance for running the walker
        """

        self.currenttask += 1

        # Give the task an identity
        task = self.wq.new_task()
        tag  = self.encode_task_tag(walker)
        task.specify_tag(tag)

        # Give the task the information it needs
        self.marshal_to_task(walker, task)
        return task

    def _try_duplicate_tasks(self):
        """
        Attempt to duplicate tasks to see if they can be run on idle workers.
        Useful in the event that some workers are faster.

        Parameters:
            None

        Returns:
            None
        """

        i = 0
        # Duplicate tasks until no more can be duplicated
        while self.wq.can_duplicate_tasks():
            i += 1
            if i > 20: break

            # Get a task that can be duplicated if one exists
            tag = self.wq.select_tag()

            # Return if no tasks can be duplicated
            if tag is None: break

            # Get walker info, create a new task for it, ans submit the task
            wid    = self.decode_from_task_tag(tag)['walkerid']
            walker = self.system.walker(wid)
            task   = self._new_task(walker)
            self.wq.submit(task)

    def _resubmit(self):
        invalid = self.system.filter_by_valid()
        while invalid != {}:
            if self.verbose:
                print("Resubmitting invalid models")
            for cid in invalid:
                valid = self.system.get_valid_walker(cid)
                if valid is None:
                    raise InvalidCellException
                for w in invalid[cid]:
                    w._start = numpy.array(valid.start())
                    task = self._new_task(w)
                    self.wq._submit(task)
            self.recv()
            invalid = self.system.filter_by_valid()

    def _recv(self):
        """
        Receive completed tasks and assign their results to the System.

        Parameters:
            None

        Returns:
            None
        """
        if self._verbose:
            print(time.asctime(), 'Receiving tasks')
        system = self.system

        # Start recording AWE statistics for the receive phase
        if self._log:
            self.stats.time_barrier('start')

        # Receive tasks until there are none left to receive
        while not self.wq.empty:
            # Get the walker results and store/save them
            walker = self.wq.recv(self.marshal_from_task, self.mark_invalid_task)
            system.set_walker(walker)
            self.logwalker(walker)

            # Attempt to duplicate any slow walkers
            self._try_duplicate_tasks()

        # Stop recording stats
        if self._log:
            self.stats.time_barrier('stop')
        self.wq.clear()
        if self._verbose:
            print(system)


    def _resample(self):
        """
        Perform resampling on the System.

        Parameters:
            None

        Returns:
            None
        """
        #print("in AWE._resample()")
        # Keep track of how long it takes to resample
        if self._log:
            self.stats.time_resample('start')

        self.system = self.resample(self.system)

        if self._log:
            self.stats.time_resample('stop')
        #print("done resampling")

    def run(self):
        """
        Run the AWE-WQ program.

        Parameters:
            None

        Returns:
            None
        """

        # If the system previously shut down, restart it at the last checkpoint
        #self.recover()

        # If this is the first run, save the initial system to file
        if self._firstrun:
            self.resample.save(self.system) # May raise NotImplementedError
            self._firstrun = False

        # Ensure that the system actually has data to run
        assert len(self.system.cells) > 0
        assert len(self.system.walkers) > 0

        # Begin recording log information
        if self._log:
            t = time.time()
            self.statslogger.update(t, 'AWE', 'start_unix_time', t)

        # Run hte specified iterations number of iterations and exit on ctrl+c
        try:
            while self.iteration < self.iterations:
                #if self._log:
                    #fd = open("memory_stats/"+float.__str__(t)+"rss_usage.csv", 'a')
                    #fd.write(int.__str__(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)+"\n")
                #print("MaxRSS Memory: %s" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
                #m = input("Press enter")

                # Update the checkpoint
                if self.iteration % self.checkpointfreq == 0:
                    if self._verbose:
                        print(time.asctime(), 'Checkpointing to', self.traxlogger.cpt_path)
                    self.checkpoint()

                # Increment the iteration
                self.iteration += 1

                # Log statistics to file
                if self._verbose:
                    print(time.asctime(), 'Iteration', self.iteration, 'with', len(self.system.walkers), 'walkers')
                if self._log:
                    runtime = stats.time.time()
                    self.statslogger.update(runtime, 'AWE', 'iteration', self.iteration)
                    self.statslogger.update(runtime, 'AWE', 'walkers', len(self.system.walkers))

                    self.stats.time_iter('start')

                # Send, receive, and resample (a.k.a. the actual work)
                self._submit()
                self._recv()     ## barrier
                self._resample()

                if self._log:
                    self.stats.time_iter('stop')


        except KeyboardInterrupt:
            pass

        # except Exception, e:
        #     print 'Failed:', e
        #     import sys
        #     sys.exit(1)


    # @typecheck(int, int)
    # @returns(str)
    def encode_task_tag(self, walker):
        """
        Encode a walker as the identifier for a task. The tag contains relevant
        information from the task to determine which walker was processed
        independent of the walker object itself (which is included as a .pkl
        file in the task).

        Parameters:
            walker - the walker to be encoded for a task

        Returns:
            The tag for a task
        """

        tag = '%(outfile)s+%(cellid)d+%(weight)f+%(walkerid)d' % {
            'outfile' : os.path.join(self.wq.tmpdir, workqueue.RESULT_NAME),
            'cellid'  : walker.assignment,
            'weight'  : walker.weight,
            'walkerid' : walker.id
        }

        return tag

    @typecheck(str)
    @returns(dict)
    def decode_from_task_tag(self, tag):
        """
        Decode walker information form the task tag.

        Parameters:
            tag - the task tag to decode

        Returns:
            A dictionary containing sufficient information to identify a walker
        """
        # print(tag)
        split = tag.split('+')
        outfile, cellid, weight, walkerid = tag.split('+')
        return {'cellid'   : int(cellid)   ,
                'weight'   : float(weight) ,
                'walkerid' : int(walkerid) ,
                'outfile'  : outfile       }



    @typecheck(int, workqueue.WQ.Task)
    def marshal_to_task(self, walker, task):
        """
        Prepare all of the necessary files to send with a task to a worker.

        Parameters:
            walker - the walker for task to run
            task   - the task to prepare for submission

        Returns:
            None
        """

        # Convert the topology to a PDB format
        top        = self.system.topology
        top.coords = walker.start
        pdbdat     = str(top)

        # Serialize the walker
        #wdat = pickle.dumps(walker).decode('raw_unicode_escape')
        #pkl_name = os.path.join(workqueue.PICKLE_BASE, "walker"+int.__str__(walker.id)+".pkl")
        #pkl_file = open(pkl_name, 'wb')
        #pickle.dump(walker, pkl_file)
        #pkl_file.close()

        # Send the the topology and walker to the worker
        # See cctools work_queue.Task for more information
        task.specify_buffer(
            pdbdat,
            #"structure.pdb",
            #sys.getsizeof(pdbdat),
            workqueue.WORKER_POSITIONS_NAME+"."+int.__str__(self.currenttask),
            cache=False
        )

        #task.specify_file(
        #    pkl_name,
        #    #sys.getsizeof(wdat),
        #    workqueue.WORKER_WALKER_NAME+"."+int.__str__(self.currenttask),
        #    0, # WORK_QUEUE_INPUT
        #    cache=False)

        self.specify_task_output_file(task)

    def specify_task_output_file(self, task):
        """
        Tell the task where to place its output.

        Parameters:
            task - the task for which to set an output dir

        Returns:
            None
        """

        output = os.path.join(self.wq.tmpdir, task.tag)
        task.specify_output_file(output, remote_name = workqueue.WORKER_RESULTS_NAME+"."+int.__str__(self.currenttask), cache=False)

    @typecheck(workqueue.WQ.Task)
    @returns(Walker)
    def marshal_from_task(self, result):
        """
        Get files returned from a task.

        Parameters:
            result - a completed task containing output

        Returns:
            The walker that was sent along with the task updated to include
            result information
        """

        # The output file is compressed, so untar it and get the relevant info
        import tarfile
        #if tarfile.is_tarfile(result.tag):
        #    print("%s is a valid tarfile" % result.tag)
        #else:
        #    print("Invalid tarfile: %s" % result.tag)
        tag_info = self.decode_from_task_tag(result.tag)
        print(tag_info)
        tar = tarfile.open(result.tag)
        try:
            #walkerstr         = tar.extractfile(workqueue.WORKER_WALKER_NAME).read()
            #print(walkerstr)

            pdbstring         = tar.extractfile(workqueue.RESULT_POSITIONS).read()
            # print(pdbstring)

            cellstring        = tar.extractfile(workqueue.RESULT_CELL).read()
            # print(cellstring)
        finally:
            tar.close()

        # Get the coordinates from the ending configuration of the walker
        pdb               = structures.PDB(pdbstring.decode("utf-8"))
        # print("Got the pdb")
        coords            = pdb.coords
        # print("Got the coordinates")

        # Get the cell id for the ending coordinates
        cellid            = int(cellstring.decode("utf-8"))
        # print("Got the cell id")

        # Load the walker object from the .pkl file
        # print(tag_info["walkerid"])
        walker            = self.system.walker(tag_info["walkerid"])#pickle.loads(walkerstr)
        # print("Got the walker from system")

        # Determine whether or not the walker changed states
        transition = walker.assignment != cellid
        # print("The transition was %s" % (transition))
        if walker is not None:
            if self._verbose:
                print(time.asctime(), 'Iteration', self.iteration, '/', self.iterations, \
                      'Walker', walker.id, \
                      'transition', walker.assignment, '->', cellid, \
                      self.wq.tasks_in_queue(), 'tasks remaining')
            #pkl_name = "walker"+int.__str__(walker.id)+".pkl"
            #os.remove(workqueue.PICKLE_BASE+pkl_name)
            # Log the walker
            if self._log:
                self.transitionslogger.update(time.time(), 'AWE', 'cell_transition',
                                      'iteration %s from %s to %s %s' % \
                                          (self.iteration, walker.assignment, cellid, transition))

            # Update the walker's state
            walker.end        = coords
            walker.assignment = cellid

        os.unlink(result.tag)
        return walker

    @typecheck(workqueue.WQ.Task)
    def mark_invalid_task(self, task):
        tag_info = self.decode_from_task_tag(task.tag)
        walker = self.system.walker(tag_info["walkerid"])
        walker.mark_invalid()





class Cell(object):

    """
    Container for the state of a Cell in the molecular state space.

    Fields:
        id    - the id of the cell
        core  - the color of the cell

    Methods:
        None
    """

    def __init__(self, cid, weight=1., core=DEFAULT_CORE, walkers=None):

        """
        Initialize a new instance of Cell.

        Parameters:
            cid     - the id of the cell
            weight  - unused; possibly legacy for MSM
            core    - the color of the cell (see resample.OneColor for usage)
            walkers - unused; possibly legacy for cells keeping a walker list
        """

        self._id      = cid
        self._core    = core


    @property
    def id(self): return self._id

    @property
    def core(self): return self._core

    @property
    def color(self, wid): return self._walkers[wid].color

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
    """
    Contains all information for a working Weighted Ensemble system.

    Fields:
        topology - the topology of the molecule represented by walkers
        cells    - the cells within the system
        walkers  - the walkers managed by the system

    Methods:
        add_cell        - add a cell to the system
        set_cell        - set the state of a cell in the system
        walker          - get a particular walker from the system
        add_walker      - add a walker to the system
        set_walker      - set the state of a walker in the system
        cell            - get a particular cell from the system
        has_cell        - determine if a cell is in the system
        filter_by_cell  - filter the walker list by walker cell id
        filter_by_color - filter the walker list by walker color
        filter_by_core  - filter the cell list by cell core
        clone           - make a copy of the system
    """

    def __init__(self, topology=None, cells=None):
        """
        Initialize a new instance of System.

        Parameters:
            topology - the topology of the molecule represented by walkers
            cells    - a dictionary of cells in the system mapped by cell id

        Returns:
            None
        """

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
        return list(self._cells.values())

    @property
    @returns(list)
    def walkers(self):
        return list(self._walkers.values())

    @property
    # @returns(np.array)
    def weights(self):
        """
        Get the weights of all walkers in the System.
        """

        return np.array(list(map(lambda w: w.weight, self.walkers)))

    @property
    @returns(set)
    def colors(self):
        """
        Get the colors of all walkers in the System.
        """

        return set([w.color for w in self.walkers])#set(map(lambda w: w.color, self.walkers))

    @typecheck(Cell)
    def add_cell(self, cell):
        """
        Add a cell to the system cell dictionary.

        Parameters:
            cell - a cell to add that has the 'id' attribute

        Returns:
            None
        """

        if cell.id in self._cells:
            raise ValueError('Duplicate cell id %d' % cell.id)
        self.set_cell(cell)

    @typecheck(Cell)
    def set_cell(self, cell):
        """
        Update the state of a cell in the System.

        Parameters:
            cell - the cell to update that has the 'id' attribute

        Returns:
            None
        """

        self._cells[cell.id] = cell

    @returns(Walker)
    def walker(self, wid):
        """
        Find a walker in the System.

        Parameters:
            wid - the id of the walker to find

        Returns:
            The walker with the supplied id or None if it is not in the System
        """

        return self._walkers[wid]

    @typecheck(Walker)
    def add_walker(self, walker):
        """
        Add a walker to the System.

        Parameters:
            walker - a walker with attribute 'id' to add

        Returns:
            None
        """

        assert walker.assignment >= 0, 'is: %s' % walker.assignment
        self.set_walker(walker)

    @typecheck(Walker)
    def set_walker(self, walker):
        """
        Update the state of a walker in the System.

        Parameters:
            walker - a walker with attribute 'id'

        Returns:
            None
        """

        self._walkers[walker.id] = walker

    @returns(Cell)
    def cell(self, i):
        """
        Get a cell from the System.

        Parameters:
            i - the id of the cell to get

        Returns:
            The cell with the supplied id or None if it is not in the System
        """

        return self._cells[i]

    @returns(bool)
    def has_cell(self, i):
        """
        Determine if a cell is in the System.

        Parameters:
            i - the id of the cell to find

        Returns:
            A Boolean value representing whether or not a call with the
            supplied id is in the System
        """

        return i in self._cells

    # @returns(System)
    def filter_by_cell(self, cell):
        """
        Filter the walker list to include only walkers in the specified cell.

        Parameters:
            cell - the cell by which to filter

        Returns:
            A new System instance containing only walkers that are in the
            supplied cell
        """

        ws     = list([w for w in self.walkers if w.assignment == cell.id])#filter(lambda w: w.assignment == cell.id, self.walkers)
        newsys = self.clone(cells={cell.id:self.cell(cell.id)})
        for w in ws: newsys.add_walker(w)
        return newsys

    # @returns(System)
    def filter_by_color(self, color):
        """
        Filter the walker list to include only walkers of a specified color.

        Parameters:
            color - the color by which to filter

        Returns:
            A new System instance containing only walkers that are of the
            supplied color
        """

        ws     = [w for w in self.walkers if w.color == color]#filter(lambda w: w.color == color, self.walkers)
        newsys = self.clone()

        for w in ws:
            newsys.add_walker(w)

            cell = self.cell(w.assignment)
            newsys.set_cell(cell)

        return newsys


    # @returns(System)
    def filter_by_core(self, core):
        """
        Filter the cell list to include only walkers of a specified core.

        Parameters:
            core - the core by which to filter

        Returns:
            A new System instance containing only cells that are of the
            supplied core
        """

        cells  = [c for c in self.cells if c.core == core]#filter(lambda c: c.core == core, self.cells)
        cs     = {}
        for c in cells: cs[c.id] = c

        newsys = self.clone(cells=cs)
        for w in self.walkers:
            if w.assignment in cs:
                newsys.add_walker(w)

        return newsys

    def filter_by_valid(self):
        invalid = {}
        for c in self.cells:
            walkers = self.filter_by_cell(c)
            if len(walkers) > 0:
                invalid[c.id] = [w for w in walkers if not w.valid]
        return invalid

    def get_valid_walker(self, cid):
        c = self.cell(cid)
        walkers = [w for w in self.filter_by_cell(c) if w.valid]
        return walkers[random.randint(0,len(walkers))] if len(walkers) > 0 else None

    def clone(self, cells=True):
        """
        Make a copy the System instance.

        Parameters:
            cells - True to copy the current cell dictionary, False to start
                    with an empty cell dictionary

        Returns:
            A new instance of System with the same topology as the caller and
            either a copy of the cell dictionary or an empty cell dictionary
            depending on the supplied flag
        """

        _cells = self._cells if cells else dict()
        return System(topology=self.topology, cells=_cells)


class SinkStates(object):
    """
    Manager for associating a color with states and vice-versa.

    Fields:
        None

    Methods:
        add    - add states by color
        color  - get the color of a cell
        states - get the set of states of a color
    """

    def __init__(self):
        """
        Initialize a new instance of SinkStates.

        Parameters:
            None

        Returns:
            None
        """

        # Note that each color has an associated ***set*** of states
        self._color_state = defaultdict(set)
        self._state_color = dict()

    def add(self, color, *states):
        """
        Add states and associate them with a color.

        Parameters:
            color  - the color to associate with the states
            states - a list of states (cell ids) to add

        Returns:
            None
        """

        for state in states:
            self._color_state[color].add(state)
            self._state_color[state] = color

    def color(self, cell):
        """
        Determine the color of a particular cell.

        Parameters:
            cell - a cell with property 'id' to find

        Returns:
            The color of the cell or the global default color if it does not
            have an associated color.
        """

        if cell.id in self._state_color:
            return self._state_color[cell.id]
        else:
            global _DEFAULT_COLOR
            return _DEFAULT_COLOR

    def states(self, color):
        """
        Get the set of states associated with a color.

        Parameters:
            color - the color to find

        Returns:
            The set of states associated with the supplied color or None if
            no states of that color exist.
        """

        return self._color_state[color]

    @property
    def ncolors(self): return len(self._color_state)
