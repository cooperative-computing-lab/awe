# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


import awe

import work_queue as WQ

import os, tarfile, tempfile, time, shutil, traceback, random
from collections import defaultdict


### A process can only support a single WorkQueue instance
_AWE_WORK_QUEUE = None


### workaround for now.
##+ These are the names of the input/output filess to be materialized on the worker

WORKER_POSITIONS_NAME     = 'structure.pdb'
WORKER_WALKER_NAME  = 'walker.pkl'
WORKER_WEIGHTS_NAME = 'weight.dat'
WORKER_COLOR_NAME   = 'color.dat'
WORKER_CELL_NAME    = 'cell.dat'
WORKER_RESULTS_NAME = 'results.tar'

RESULT_POSITIONS    = 'structure2.pdb'
RESULT_WEIGHTS      = 'weight.dat'
RESULT_COLOR        = 'color.dat'
RESULT_CELL         = 'cell2.dat'
RESULT_NAME         = 'results.tar'


class WorkQueueException       (Exception): pass
class WorkQueueWorkerException (Exception): pass

class WQFile(object):

    @awe.typecheck(str, base=bool, cached=bool)
    def __init__(self, masterpath, base=True, cached=True):
        self._masterpath = masterpath
        self._base       = base
        self._cached     = cached

    @property
    def masterpath(self):
        return self._masterpath

    @property
    def remotepath(self):
        if self.isbase:
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
        if '$' not in self.masterpath and not os.path.exists(self.masterpath):
            raise IOError, 'Cannot find file to send to worker: %s' % self.masterpath
        task.specify_file(self.masterpath, remote_name=self.remotepath, cache=self.cached)

    def __str__(self):
        return 'WQFile: masterpath=%s remotepath=%s cached=%s' % (self.masterpath, self.remotepath, self.cached)

    def __repr__(self):
        return 'WQFile(%r, base=%r, cached=%r' % (self._masterpath, self._base, self._cached)



class Config(object):
    """
    Class for configuring a WorkQueue instance
    """

    def __init__(self):

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
        self.monitor         = False	
        self.summaryfile     = ''
	self.capacity        = False

        self._executable = None
        self._cache = set()

    executable = property(lambda self: self._executable)
    getcache   = property(lambda self: self._cache)

    def execute(self, path):
        f = WQFile(path)
        self._executable = f
        self.cache(f.masterpath)

    def cache(self, *files, **kws):
        base = kws.get('base', True)
        for path in files:
            wqf = WQFile(path, base=base, cached=True)
            self._cache.add(wqf)

    def _mk_wq(self):
        global _AWE_WORK_QUEUE
        if _AWE_WORK_QUEUE is not None:
            ### warn
            awe.log('WARNING: using previously created WorkQueue instance')
        else:
            if self.debug:
                WQ.set_debug_flag(self.debug)
                if self.wq_logfile:
                     awe.util.makedirs_parent(self.wq_logfile)
                     WQ.cctools_debug_config_file(self.wq_logfile)
                     WQ.cctools_debug_config_file_size(0) 
            if self.name:
                self.catalog = True
            wq = WQ.WorkQueue(name      = self.name,
                              port      = self.port,
                              shutdown  = self.shutdown,
                              catalog   = self.catalog,
                              exclusive = self.exclusive)
            wq.specify_algorithm(self.schedule)
            if self.monitor: 
                wq.enable_monitoring(self.summaryfile)

	    if self.capacity:
		wq.estimate_capacity()
 
            awe.log('Running on port %d...' % wq.port)
            if wq.name:
                awe.log('Using project name %s' % wq.name)
            if self.debug and self.wq_logfile:
                awe.log('Logging WorkQueue to %s' % self.wq_logfile)

            typ = type(self.fastabort)
            if typ is float or typ is int:
                wq.activate_fast_abort(self.fastabort)

            _AWE_WORK_QUEUE = wq

        awe.util.makedirs_parent(self.wqstats_logfile)
        _AWE_WORK_QUEUE.specify_log(self.wqstats_logfile)
        return _AWE_WORK_QUEUE


class TagSet(object):
    def __init__(self, maxreps=5):
        self._tags = defaultdict(set)
        self._maxreps = maxreps

    def can_duplicate(self):
        valid = filter(lambda k: k < self._maxreps, self._tags.iterkeys())
        return len(valid) > 0

    def clear(self):
        self._tags.clear()

    def clean(self):
        for k in self._tags.keys():
            if len(self._tags[k]) < 1:
                del self._tags[k]

    def _find_tag_group(self, tag):
        for group, tags in self._tags.iteritems():
            if tag in tags:
                return group
        return None

    def add(self, tag, startcount=0):
        key = self._find_tag_group(tag)

        ### add the tag to the appropriate group, removing it from previous one
        if key is None:
            self._tags[startcount].add(tag)
        else:
            self._tags[key+1].add(tag)
            self._tags[key  ].discard(tag)

        ### delete the group if it became empty
        if key is not None and len(self._tags[key]) == 0:
            del self._tags[key]


    def select(self):
        if len(self) > 0:
            count  = 1
            minkey = min(self._tags.keys())
            assert len(self._tags[minkey]) > 0, str(minkey) + ', ' + str(self._tags[minkey])
            return random.sample(self._tags[minkey], count)[0]
        else:
            return None

    def discard(self, tag, key=None):
        key = key or self._find_tag_group(tag)
        if key is not None:
            self._tags[key].discard(tag)

    def __len__(self):
        return reduce(lambda s, k: s + len(self._tags[k]), self._tags.iterkeys(), 0 )

    def __str__(self):
        d = dict([(k,len(s)) for k,s in self._tags.iteritems()])
        return '<TagSet(maxreps=%s): %s>' % (self._maxreps, d)



class WorkQueue(object):

    # @awe.typecheck(Config)
    def __init__(self, cfg, statslogger=None, taskoutputlogger=None):

        self.cfg    = cfg
        self.wq     = self.cfg._mk_wq()
        self._tagset = TagSet(maxreps=self.cfg.maxreps)

        self.stats  = awe.stats.WQStats(logger=statslogger)

        self.tmpdir = tempfile.mkdtemp(prefix='awe-tmp.')

        self.restarts = dict()

        self.statslogger      = statslogger      or awe.stats.StatsLogger(buffersize=42)
        self.taskoutputlogger = taskoutputlogger or awe.stats.StatsLogger(path='debug/task_output.log.gz', buffersize=42)

    @property
    def empty(self):
        return self.wq.empty()

    def __getstate__(self):
        """
        SwigPyObjects cannot be pickles, so remove the underlying WorkQueue object
        """
        self.statslogger.close()
        self.taskoutputlogger.close()
        odict = self.__dict__.copy()
        del odict['wq']
        return odict

    def __setstate__(self, odict):
        """
        Since SwigPyObjects are not pickleable, we just recreate the WorkQueue object from the configuration
        """
        self.__dict__.update(odict)
        self.statslogger.open()
        self.taskoutputlogger.open()


    def __del__(self):
        import shutil
        shutil.rmtree(self.tmpdir)


    @awe.typecheck(WQ.Task)
    def update_task_stats(self, task):
        self.stats.task(task)

    def new_task(self):
        cmd = self.cfg.executable.remotepath
        task = WQ.Task('./' + cmd)

        ### executable
        self.cfg.executable.add_to_task(task)

        ### cached files
        for wqf in self.cfg.getcache:
            wqf.add_to_task(task)

        return task

    @awe.typecheck(WQ.Task)
    def submit(self, task):
        self._tagset.add(task.tag)
        return self.wq.submit(task)

    @awe.typecheck(WQ.Task)
    def restart(self, task):
        if task.tag not in self.restarts:
            self.restarts[task.tag] = 0

        if self.restarts[task.tag] < self.cfg.restarts:
            print time.asctime(), 'task failed with', task.return_status, \
                'result is', task.result, \
                'restarting', task.tag, \
                '#%d' % (self.restarts[task.tag] + 1)
            self.submit(task)
            self.restarts[task.tag] += 1
            return True
        else:
            return False

    def wait(self, *args, **kws):
        return self.wq.wait(*args, **kws)

    def taskoutput(self, task):
        output = task.output or ''
        output = ('\n' + output).split('\n')
        output = '\n\t'.join(output)
        return output

    def add_tag(self, tagtext):
        self._tagset.add(tagtext)

    def discard_tag(self, tagtext):
        self._tagset.discard(tagtext)

    def cancel_tag(self, tagtext):
        while self.wq.cancel_by_tasktag(tagtext):
            pass

    def select_tag(self):
        self._tagset.clean()
        return self._tagset.select()

    def clear_tags(self):
        self._tagset.clear()

    def tasks_in_queue(self):
        return self.wq.stats.tasks_running + self.wq.stats.tasks_waiting

    def active_workers(self):
        return self.wq.stats.workers_busy + self.wq.stats.workers_ready 

    def can_duplicate_tasks(self):
        return  self.tasks_in_queue() < self.active_workers() \
            and self._tagset.can_duplicate()


    def recv(self, marshall):

        # print time.asctime(), 'waiting for task'
        while True:

            task = self.wait(self.cfg.waittime)

            if task:
                self.update_task_stats(task)
                # print time.asctime(), 'Received result. %d tasks remaining in iteration.' % self.tasks_in_queue()

                self.taskoutputlogger.output("<====== WQ: START task %s output ======>\n" % task.tag)
                self.taskoutputlogger.output(task.output)
                self.taskoutputlogger.output("<====== WQ: END task %s output ======>\n"   % task.tag)

            if task and task.result == 0:

                if not task.return_status == 0 and not self.restart(task):
                    raise WorkQueueWorkerException, \
                        self.taskoutput(task) + '\n\nTask %s failed with %d' % (task.tag, task.return_status)


                try:
                    result = marshall(task)
                except Exception, ex:

                    ### sometimes a task fails, but still returns.
                    ##+ attempt to restart these
                    if not self.restart(task):
                        raise WorkQueueException, \
                            self.taskoutput(task) + '\n\nMaster failed: could not load resultfile:\n %s: %s\n\n%s' % \
                            (ex.__class__.__name__, ex, traceback.format_exc())
                    else:
                        continue

                self.cancel_tag(task.tag)
                self.discard_tag(task.tag)
                return result

            elif task and not task.result == 0:
                if not self.restart(task):
                    raise WorkQueueException, 'Task exceeded maximum number of resubmissions for %s\n\n%s' % \
                        (task.tag, self.taskoutput(task))

            else: continue
