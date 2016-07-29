# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


import awe

import work_queue as WQ

import os, tarfile, tempfile, time, shutil, traceback, random, re
from collections import defaultdict
from functools import reduce

### A process can only support a single WorkQueue instance
_AWE_WORK_QUEUE = None


### workaround for now.
##+ These are the names of the input/output filess to be
##+ materialized on the worker

# These seem to have been deprecated in production.
# The workers use a different set of files.

PICKLE_BASE = 'awe-instance-data/pickle/'

WORKER_POSITIONS_NAME = 'structure.pdb' # The PDB used to generate a trajectory
WORKER_WALKER_NAME    = 'walker.pkl'    # The walker object for the task
WORKER_WEIGHTS_NAME   = 'weight.dat'    # The weight for each walker
WORKER_COLOR_NAME     = 'color.dat'     # The color (state?) of the walker
WORKER_CELL_NAME      = 'cell.dat'      # The list of exemplar configurations
WORKER_RESULTS_NAME   = 'results.tar'   # The file to return to the master

RESULT_POSITIONS = 'structure2.pdb' # The PDB of the final trajectory frame
RESULT_WEIGHTS   = 'weight.dat'     # The weight of the resulting configuration
RESULT_COLOR     = 'color.dat'      # The color (state) of the walker
RESULT_CELL      = 'cell2.dat'      # The updated cells information
RESULT_NAME      = 'results.tar'    # The file to return to the master


class WorkQueueException       (Exception): pass
class WorkQueueWorkerException (Exception): pass

class WQFile(object):
    """
    Manages filepath information and options for temp directories containing
    files used by cctools WorkQueue.Task objects.

    Fields:
        _masterpath - the path to the master set of files
        _base       - flag for looking in the current directory
        _cached     - flag for determining if the path points to cached
                      task files
        _remotepath - the path to the master set of files if not local

    Methods:
        add_to_task - tell a task where to find the files it needs
    """

    @awe.typecheck(str, base=bool, cached=bool, remotepath=str)
    def __init__(self, masterpath, base=True, cached=True, remotepath=None):
        """
        Initialize a new instance of WQFile.

        Properties:
            masterpath - the path to the master set of task files
            base       - whether or not the current directory is the directory
                         to look in
            cached     - whether or not the path points to files cached for a
                         WorkQueue Task
            remotepath - the remote directory to look in for task files
        """

        self._masterpath = masterpath
        self._base       = base
        self._cached     = cached
        self._remotepath = remotepath

    @property
    def masterpath(self):
        return self._masterpath

    @property
    def remotepath(self):
        if self._remotepath:
            return self._remotepath
        elif self.isbase:
            # This returns only a filename, look in the current directory
            return os.path.basename(self.masterpath)
        else:
            return self.masterpath

    @property
    def isbase(self):
        return self._base

    @property
    def cached(self):
        return self._cached

    def add_to_task(self, task):
        """
        Tell a task where to find the files it needs to do its work.

        Parameters:
            task - a cctools WorkQueue.Task object

        Returns:
            None
        """

        if '$' not in self.masterpath and not os.path.exists(self.masterpath):
            raise IOError('Cannot find file to send to worker: %s' % self.masterpath)
        task.specify_file(self.masterpath, remote_name=self.remotepath, cache=self.cached)

    def __str__(self):
        return 'WQFile: masterpath=%s remotepath=%s cached=%s' % (self.masterpath, self.remotepath, self.cached)

    def __repr__(self):
        return 'WQFile(%r, base=%r, cached=%r' % (self._masterpath, self._base, self._cached)



class Config(object):
    """
    Configuration options for a WorkQueue instance.

    Fields:
        name            - the name of this run of AWE-WQ
        port            - the port on which the WorkQueue master listens
        schedule        - the task scheduling algorithm to use
        exclusive       - whether or not the WorkQueue instance is singleton
        catalog         -
        debug           - which information to include in logs
        shutdown        -
        fastabort       -
        restarts        - the number of times to restart a failed task
        maxreps         -
        waittime        -
        wq_logfile      - the log file for WorkQueue debug information
        wqstats_logfile - the log file for WorkQueue statistical information
        monitor         -
        summaryfile     - the log file for WorkQueue run basic information
        capacity        -
        _executable     - the script or executable that the task should run
        _cache          -

    Methods:
        execute - add the main program to run to the task cached file list
        cache   - add files to the task cached file list
        _mk_wq  - create a WorkQueue instance or get a reference to the
                  singleton
    """

    def __init__(self):
        """
        Initialize a new instance of Config. Uses default options that are then
        reset by the aweclasses.AWE class if necessary.

        Parameters:
            None

        Returns:
            None
        """

        self.name            = ''
        self.port            = WQ.WORK_QUEUE_DEFAULT_PORT
        self.schedule        = WQ.WORK_QUEUE_SCHEDULE_TIME
        self.exclusive       = True
        self.catalog         = False
        self.debug           = ''
        self.shutdown        = False
        self.fastabort       = 3
        self.restarts        = 95 # until restarts are handled on a per-iteration basis
        self.maxreps         = 9
        self.waittime        = 10 # in seconds
        self.wq_logfile      = 'debug/wq.log'
        self.wqstats_logfile = 'debug/wq-stats.log'
        self.wqtransactions_logfile = 'debug/wq-transactions.log'
        self.monitor         = True	
        self.summaryfile     = 'debug/monitor-summaries'
        self.capacity        = False
        self.task_config     = {
                                   "cores": 1,
                                   "memory": 250,
                                   "disk": 250,
                               }
        self._executable = None
        self._cache = set()

    executable = property(lambda self: self._executable)
    getcache   = property(lambda self: self._cache)

    def execute(self, path):
        """
        Tag a file as the intended program for the task to execute and add it
        to the cache.

        Parameters:
            path - filepath of the program to execute

        Returns:
            None
        """

        f = WQFile(path)
        self._executable = f
        self._cache.add(f)

    def cache(self, *files, **kws):
        """
        Add a list of files to the WorkQueue task cache. Cached files are saved
        between tasks run by a single worker.

        Parameters:
            files - a tuple of filepaths
            kws   - a dictionary of keyword-mapped arguments

        Keyword Arguments:
            base       - a boolean value representing whether to use the
                         current directory
            remotepath - the remote filepath to look in for files

        Returns:
            None
        """

        base = kws.get('base', True)
        remotepath = kws.get('remotepath', '')
        for path in files:
            wqf = WQFile(path, base=base, cached=True, remotepath=remotepath)
            self._cache.add(wqf)

    def _mk_wq(self):
        """
        Only one instance of WorkQueue should be run per process. This grants
        access to the WorkQueue singleton or else creates a new WorkQueue
        instance. This also ensures that the cctools WorkQueue object can
        handle more workers.

        Parameters:
            None

        Returns:
            The cctools WorkQueue singleton object
        """

        global _AWE_WORK_QUEUE
        if _AWE_WORK_QUEUE is not None:
            ### warn
            awe.log('WARNING: using previously created WorkQueue instance')
        else:
            if self.debug:
                # Set up debugging parameters for the cctools WorkQueue object.
                # It has inbuilt debugging capabilities.
                WQ.set_debug_flag(self.debug)

                if self.wq_logfile:
                     awe.util.makedirs_parent(self.wq_logfile)
                     WQ.cctools_debug_config_file(self.wq_logfile)
                     WQ.cctools_debug_config_file_size(0)

            if self.name:
                self.catalog = True

            # Create the cctools WorkQueue object
            wq = WQ.WorkQueue(name      = self.name,
                              port      = self.port,
                              shutdown  = self.shutdown,
                              catalog   = self.catalog,
                              exclusive = self.exclusive)

            # Specify the task scheduling algorithm
            wq.specify_algorithm(self.schedule)

            # Turn cctools WorkQueue object status monitoring on or off
            if self.monitor:
                wq.enable_monitoring(self.summaryfile)
                wq.specify_transactions_log(self.wqtransactions_logfile)

            if self.capacity:
                # Determine the number of workers the WorkQueue object can handle
    	        wq.estimate_capacity()

            # Display information about this run of AWE-WQ
            awe.log('Running on port %d...' % wq.port)
            if wq.name:
                awe.log('Using project name %s' % wq.name)
            if self.debug and self.wq_logfile:
                awe.log('Logging WorkQueue to %s' % self.wq_logfile)

            # Set up fast abort procedures
            typ = type(self.fastabort)
            if typ is float or typ is int:
                wq.activate_fast_abort(self.fastabort)

            # Ensure that the singleton is set to the new instance
            _AWE_WORK_QUEUE = wq

        # Ensure that the singleton is logging to the correct files
        awe.util.makedirs_parent(self.wqstats_logfile)
        _AWE_WORK_QUEUE.specify_log(self.wqstats_logfile)

        # Return a reference to the singleton
        return _AWE_WORK_QUEUE


class TagSet(object):
    """
    Manages tags for identifying submitted tasks. Tags are stored in sets
    indexed in a dictionary by the total number of duplicates existing for
    that particular tag (e.g. if tag in self._tags[3], there exist three
    duplicates of that tag).

    Fields:
        _tags    - the dictionary of sets of tags
        _maxreps - the maximum number of times a tag may be duplicated

    Methods:
        can_duplicate   - determines whether any task can be duplicated
        clear           - clear the tag set dictionary
        clean           - remove a tag set from the dictionary if it contains
                          no tags
        _find_tag_group - find the set containing a specific tag in the
                          dictionary
        add             - add a tag (or duplicate) to the dictionary
        select          - get a random tag from the tag set with the least
                          duplicates
        discard         - remove a tag from the dictionary
    """

    def __init__(self, maxreps=5):
        """
        Initialize a new instance of TagSet.

        Parameters:
            maxreps - the maximum number of times a tag may be duplicated

        Returns:
            None
        """

        # Return an empty set if a key does not exist (see Python docs)
        self._tags = defaultdict(set)
        self._maxreps = maxreps

    def can_duplicate(self):
        """
        Determine whether any tags can can be duplicated.

        Parameters:
            None

        Return:
            A Boolean value representing whether any tag can be duplicated
        """

        #valid = filter(lambda k: k < self._maxreps, self._tags.iterkeys())
        valid = [k for k in iter(self._tags.keys()) if k < self._maxreps]
        return len(valid) > 0

    def clear(self):
        """
        Remove all tag sets from the dictionary.

        Parameters:
            None

        Returns:
            None
        """

        self._tags.clear()

    def clean(self):
        """
        Remove any empty tag sets from the dictionary.

        Parameters:
            None

        Returns:
            None
        """

        for k in list(self._tags.keys()):
            if len(self._tags[k]) < 1:
                del self._tags[k]

    def _find_tag_group(self, tag):
        """
        Find the set that the supplied tag resides in.

        Parameters:
            tag - the tag to search the dictionary for

        Returns:
            The set containing the tag or None if the tag is not in the
            dictionary
        """

        for group, tags in self._tags.items():
            if tag in tags:
                return group
        return None

    def add(self, tag, startcount=0):
        """
        Adds a tag to the dictionary in the set representing the number of
        duplicates of that tag.

        Parameters:
            tag        - the tag to add to the dictionary
            startcount - the set index to begin at

        Returns:
            None
        """

        key = self._find_tag_group(tag)

        ### add the tag to the appropriate group, removing it from previous one
        if key is None:
            # Note: this would only be dangerous if using a list
            # Dictionary keys need not be positive (or numeric), however given
            # the definition of the contents of the TagSet, this method leaves
            # the option available to have negative duplicates of a tag, etc.,
            # even though that does not make sense and would not inherently
            # cause problems.
            self._tags[startcount].add(tag)
        else:
            # Add a duplicate of the tag and remove it from its original set
            self._tags[key+1].add(tag)
            self._tags[key  ].discard(tag)

        ### delete the group if it became empty
        if key is not None and len(self._tags[key]) == 0:
            # Clean up any sets left empty
            del self._tags[key]


    def select(self):
        """
        Gets a random tag from the set of tags representing the lowest number
        of duplicates.

        Parameters:
            None

        Returns:
            A random tag from the lowest-indexed tag set or None of the
            dictionary is empty.
        """

        if len(self) > 0:
            count  = 1

            # Keys are ordinal (likely integers), so min(keys) is the tag set
            # representing tags with the least number of duplicates.
            minkey = min(self._tags.keys())

            assert len(self._tags[minkey]) > 0, str(minkey) + ', ' + str(self._tags[minkey])

            return random.sample(self._tags[minkey], count)[0]

        else:
            return None

    def discard(self, tag, key=None):
        """
        Remove a tag from the dictionary completely (as opposed to reducing the
        number of duplicates).

        Parameters:
            tag - the tag to remove
            key - the dictionary key to check, None if TagSet should determine
                  which set contains the tag

        Returns:
            None
        """

        key = key or self._find_tag_group(tag)
        if key is not None:
            # This is safe in the event that the tag is not at the
            # user-supplied key. Discard will do nothing if the element to get
            # rid of is not in the set.
            self._tags[key].discard(tag)

    def __len__(self):
        """
        Override of the object __len__ function.

        Parameter:
            None

        Returns:
            The total number of tags in the dictionary
        """

        return reduce(lambda s, k: s + len(self._tags[k]), iter(self._tags.keys()), 0 )

    def __str__(self):
        d = dict([(k,len(s)) for k,s in self._tags.items()])
        return '<TagSet(maxreps=%s): %s>' % (self._maxreps, d)



class WorkQueue(object):
    """
    An interface to the cctools work_queue module.

    Fields:
        cfg              - configuration settings container
        wq               - cctools work_queue.WorkQueue object
        _tagset          - the TagSet object containing task tags
        stats            - statistical unit from stats module
        tmpdir           - temp directory where task information is stored
        restarts         - a dictionary of tasks that have been restarted
        statslogger      - logging unit for global AWE-WQ statistics
        taskoutputlogger - logging unit for individual task output and stats

    Methods:
        update_task_stats   -
        new_task            -
        submit              -
        restart             -
        wait                -
        taskoutput          -
        add_tag             -
        discard_tag         -
        cancel_tag          -
        select_tag          -
        clear_tags          -
        clear               -
        tasks_in_queue      -
        active_workers      -
        can_duplicate_tasks -
        recv                -
    """

    # @awe.typecheck(Config)
    def __init__(self, cfg, statslogger=None, taskoutputlogger=None, verbose=False, log_it=False):
        """
        Initialize a new instance of WorkQueue for managing the cctools
        work_queue.WorkQueue object with the necessary parameters for running
        AWE-WQ.

        Parameters:
            cfg              - an initialized workqueue.Config object to use
            statslogger      - a logger object for logging statistics about
                               WorkQueue
            taskoutputlogger - a logger to use for logging the output of
                               individual tasks

        Returns:
            None
        """

        # Configure and gain access to the cctools WorkQueue singleton
        self.cfg    = cfg
        self.wq     = self.cfg._mk_wq()

        # Create the task tag set for managing duplicate tasks
        self._tagset = TagSet(maxreps=self.cfg.maxreps)

        # Start collecting statistics on the WorkQueue object
        self.stats  = awe.stats.WQStats(logger=statslogger)

        # Create the temporary file where the task cache is stored
        self.tmpdir = str(bytes(tempfile.mkdtemp(prefix='awe-tmp.'), 'ASCII'), 'ASCII')
        #print("Temporary directory is %s" % self.tmpdir)
        # Create a dictionary to keep track of the number of times a
        # particular task has been restarted
        self.restarts = dict()

        # Start logging statistics about and output from the tasks
        self.statslogger      = statslogger      or awe.stats.StatsLogger(buffersize=42)
        self.taskoutputlogger = taskoutputlogger or awe.stats.StatsLogger(path='debug/task_output.log.gz', buffersize=42)
        self._log = log_it
        self.verbose = verbose

    @property
    def empty(self):
        return self.wq.empty()

    def __getstate__(self):
        """
        SwigPyObjects cannot be pickled, so remove the underlying WorkQueue object
        See pickle documentation for more information on __getstate__
        """
        self.statslogger.close()
        self.taskoutputlogger.close()
        odict = self.__dict__.copy()
        del odict['wq']
        return odict

    def __setstate__(self, odict):
        """
        Since SwigPyObjects are not pickleable, we just recreate the WorkQueue object from the configuration
        See pickle documentation for more information on __setstate__
        """
        self.__dict__.update(odict)
        self.statslogger.open()
        self.taskoutputlogger.open()


    def __del__(self):
        """
        Standard Python delete function. It does occasionally cause problems
        if the garbage collector does not want to cooperate. If so, add tmpdir
        to the WorkQueue instance in the AWE object in aweclasses after
        initializing the WorkQueue instance. This seems to overcome the issue.

        Parameters:
            None

        Returns:
            None
        """

        import shutil
        shutil.rmtree(self.tmpdir)


    @awe.typecheck(WQ.Task)
    def update_task_stats(self, task):
        """
        Add information about a task to the WorkQueue object logging utility.
        See stats.WQStats for more information on the default logger.

        Parameters:
            task - a WorkQueue task

        Returns:
            None
        """
        if self._log:
            self.stats.task(task)

    def new_task(self):
        """
        Generate a new task object and assign it a program to run. Ensures each
        task has the correct set of supporting files. See WorkQueue.Config for
        information on task files.

        Parameters:
            None

        Returns:
            A new cctools WorkQueue.Task instance
        """

        cmd = self.cfg.executable.remotepath
        task = WQ.Task('./' + cmd)
        task.specify_cores(self.cfg.task_config["cores"])
        task.specify_memory(self.cfg.task_config["memory"])
        task.specify_disk(self.cfg.task_config["disk"])
        ### executable
        self.cfg.executable.add_to_task(task)

        ### cached files
        for wqf in self.cfg.getcache:
            wqf.add_to_task(task)

        return task

    @awe.typecheck(WQ.Task)
    def submit(self, task):
        """
        Submit a task to a worker. This adds it to the TagSet to keep track of
        duplicate tasks.

        Parameters:
            task - a configured cctools work_queue.Task object

        Returns:
            ???
        """

        self._tagset.add(task.tag)
        return self.wq.submit(task)

    @awe.typecheck(WQ.Task)
    def restart(self, task):
        """
        Restart a task if something has gone wrong and record the number of
        times that task has been restarted.

        Parameters:
            task - a configured cctools work_queue.Task object

        Return:
            A Boolean value representing whether or not the task was restarted
        """

        # Add to the dictionary of restarted tasks if first restart
        if task.tag not in self.restarts:
            self.restarts[task.tag] = 0

        # Record the restart and notify the user, then restart if the task has
        # not exceeded the maximum number of restarts
        if self.restarts[task.tag] < self.cfg.restarts:
            print(time.asctime(), 'task failed with', task.return_status, \
                'result is', task.result, \
                'restarting', task.tag, \
                '#%d' % (self.restarts[task.tag] + 1))
            self.submit(task)
            self.restarts[task.tag] += 1
            return True

        # Otherwise, do not restart the task
        else:
            return False

    def wait(self, *args, **kws):
        """
        Set the cctools work_queue.WorkQueue instance to idle state (usually
        if no workers are available).

        Parameters: (see the cctools module work_queue.py for more information)
            args - arguments for wq.wait
            kws  - keyword arguments for wq.wait

        Returns:
            ???
        """

        return self.wq.wait(*args, **kws)

    def taskoutput(self, task):
        """
        Concatenate all output from a task into a single string.

        Parameters:
            task - a configured cctools work_queue.Task object

        Returns:
            A string containing all output from the task
        """

        output = task.output or ''
        output = ('\n' + output).split('\n')
        output = '\n\t'.join(output)
        return output

    def add_tag(self, tagtext):
        """
        Add a task tag to the internal tag set dictionary.

        Parameters:
            tagtext - a string identifying the tag

        Returns:
            None
        """

        self._tagset.add(tagtext)

    def discard_tag(self, tagtext):
        """
        Discard a tag from the internal tag set dictionary.

        Parameters:
            tagtext - string identifying the tag to be removed

        Returns:
            None
        """

        self._tagset.discard(tagtext)

    def cancel_tag(self, tagtext):
        """
        Cancel all tasks associated with the supplied tag.

        Parameters:
            tagtext - a string identifying the tasks to cancel

        Returns:
            None
        """

        while self.wq.cancel_by_tasktag(tagtext):
            pass

    def select_tag(self):
        """
        Get a tag from the internal tag set dictionary. See TagSet.select for
        details on the selectoin process.

        Parameters:
            None

        Returns:
            A tag from the internal tag set dictionary
        """

        self._tagset.clean()
        return self._tagset.select()

    def clear_tags(self):
        """
        Clear the internal tag set dictionary (i.e. remove all tags).

        Parameters:
            None

        Returns:
            None
        """

        self._tagset.clear()

    def clear(self):
        """
        Remove all tags from the internal tag set dictionary and all tasks from
        the cctools work_queue.WorkQueue object instance.
        """

        self.clear_tags()

        # force clearing to allow GC, otherwise linear memory growth
        self.wq._task_table.clear()

    def tasks_in_queue(self):
        """
        The number of tasks currently running and waiting in the queue.

        Parameters:
            None

        Returns:
            An integer representing the total number of enqueued tasks
        """

        return self.wq.stats.tasks_running + self.wq.stats.tasks_waiting

    def active_workers(self):
        """
        The number of workers currently working or ready to work.

        Parameters:
            None

        Returns:
            An integer representing the number of workers not idle
        """

        return self.wq.stats.workers_busy + self.wq.stats.workers_ready

    def can_duplicate_tasks(self):
        """
        Determine whether any task can be duplicated and sent to an idle
        worker.

        Parameters:
            None

        Returns:
            A Boolean value representing whether any task can be duplicated and
            sent to an idle worker.
        """

        return  self.tasks_in_queue() < self.active_workers() \
            and self._tagset.can_duplicate()


    def recv(self, marshall, mark_invalid):
        """
        Deal with tasks as they return. Handle successes, errors, and restarts.
        Runs forever.

        Parameters:
            marshall - function describing how to send tasks to workers

        Returns:
            None
        """

        # print time.asctime(), 'waiting for task'
        while True:
            #print(self.tmpdir)
            # Set the cctools WorkQueue object to idle until a task arrives
            task = self.wait(timeout=self.cfg.waittime)

            if task:
                #print("Success!")
                # Record the task output whether it succeeded or failed
                self.update_task_stats(task)
                # print time.asctime(), 'Received result. %d tasks remaining in iteration.' % self.tasks_in_queue()
                if self._log:
                    self.taskoutputlogger.output("<====== WQ: START task %s output ======>\n" % task.tag)
                    self.taskoutputlogger.output(task.output)
                    self.taskoutputlogger.output("<====== WQ: END task %s output ======>\n"   % task.tag)
                #print(task.output)
                #m = input("press enter")
            if task and task.result == 0:
                print("Task returned with result %s and value %s" % (task.result, task.return_status) )

                # Deal with tasks in which an error occurred
                if not task.return_status == 0 and not self.restart(task):
                    raise WorkQueueWorkerException(self.taskoutput(task) + '\n\nTask %s failed with %d' % (task.tag, task.return_status))


                try:
                    result = marshall(task)
                except Exception as ex:
                    print("In the exception")
                    ### sometimes a task fails, but still returns.
                    ##+ attempt to restart these
                    if not self.restart(task):
                        raise WorkQueueException(self.taskoutput(task) + '\n\nMaster failed: could not load resultfile:\n %s: %s\n\n%s' % \
                            (ex.__class__.__name__, ex, traceback.format_exc()))
                    else:
                        continue

                # Kill the task if it cannot be restarted
                self.cancel_tag(task.tag)
                self.discard_tag(task.tag)
                return result

            elif task and not task.result == 0:
                # Check the task output for a bad model
                if task.output.find('Exception: Particle coordinate is NaN'):
                    if self.verbose:
                        print(''.join([task.tag, ': Particle coordinate is NaN']))
                    mark_invalid(task)
                    continue
                # Kill the task if it cannot be restarted.
                if not self.restart(task):
                    raise WorkQueueException('Task exceeded maximum number of resubmissions for %s\n\n%s' % \
                        (task.tag, self.taskoutput(task)))

            else: continue
