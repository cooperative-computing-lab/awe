# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


import workqueue
from util import typecheck

import numpy as np

import time as systime
import bz2


class Timer(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.t0 = 0
        self.t1 = float('inf')

    def start(self):
        self.t0 = systime.time()

    def stop(self):
        self.t1 = systime.time()

    def isrunning(self):
        return self.t0 > 0

    def elapsed(self, current=True, units='s'):
        mults      = {'s' : 1,
                      'm' : 60,
                      'h' : 60*60,
                      'd' : 60*60*24 }
        multiplier = mults[units]

        ### get the current time or from when 'Timer.stop()' was called
        if current or self.t1 == float('inf'):
            t1 = systime.time()
        else:
            t1 = self.t1

        diff       = t1 - self.t0
        return multiplier * diff


### use a single global timer for the system
_TIMER = Timer()

### simulate built-in *time* module
class time:

    @staticmethod
    def start():
        global _TIMER
        _TIMER.start()

    @staticmethod
    def time():
        global _TIMER
        t = _TIMER.elapsed(units='s')
        return t


    @staticmethod
    def timer():
        global _TIMER
        return _TIMER


class ExtendableArray(object):
    """
    A numpy.array that can be efficiently appended to
    """

    def __init__(self, typ=np.zeros, size=500, factor=2):

        self._type    = typ
        self._size0   = size
        self._count   = 0
        self._factor  = factor
        self._vals    = self._initialize()

    def _initialize(self):
        return self._type(self._size0)

    def get(self):
        return self._vals[:self._count]

    def _realloc(self):
        """
        Reallocates space if the current size and the underlying size are equal
        """

        if self._count == len(self._vals):

            x     = len(self._vals)
            alloc = x * self._factor
            self._size0 = alloc
            vals2 = self._initialize()
            vals2[:self._count] = self._vals[:self._count]
            self._vals = vals2


    def append(self, *values):

        self._realloc()
        i = self._count
        j = i + len(values)
        self._vals[i:j] = np.array(values)
        self._count += len(values)

    def __getitem__(self, i):
        if i < 0: k = self._count + i
        else:     k = i

        assert k <= self._count
        return self._vals[k]

    def __setitem__(self, i, v):
        if i < 0: k = self._count + i
        else:     k = i

        assert k <= self._count
        self._vals[k] = v

    def __len__(self):
        return self._count

    def __contains__(self, v):
        return v in self._vals[:self._count]

    def __iter__(self):
        return iter(self._vals[:self._count])

    def __str__(self):
        return str(self._vals[:self._count])

    def __repr__(self):
        return repr(self._vals[:self._count])
            

class Statistics(object):

    """
    Keep track of running statistics as the program runs
    """

    def __init__(self):

        self._num    = 0
        self._mean   = 0.
        self._m2     = 0.

        self._values = ExtendableArray() # keep track of the actual values for plots

    num  = property(lambda self: self._num)
    mean = property(lambda self: self._mean)
    var  = property(lambda self: self._m2 / float(self._num))
    std  = property(lambda self: math.sqrt(self.var))

    @property
    def values(self):
        return self._values.get()

    def update(self, *values):

        self._values.append(*values)

        for t in values:
            self._num += 1
            delta      = t - self._mean
            self._mean = self._mean + delta / float(self._num)
            self._m2   = self._m2 + delta * (t - self._mean)


    def __str__(self):
        return 'N = %s mean = %s std = %s' % (self.num, self.mean, self.std)



class WQStats(object):

    """
    Keep track of the WQ statistics
    """

    def __init__(self, logger=None):

        self.logger = logger or StatsLogger()

        self._task_times             = ExtendableArray()  # keep track of the times values are added
        self._wq_times               = ExtendableArray()

        ### task stats
        self.computation_time        = Statistics()
        self.total_bytes_transferred = Statistics()
        self.total_transfer_time     = Statistics()
        self.task_life_time          = Statistics()       # WQ Task.finish_time - Task.submit_time



    @typecheck(workqueue.WQ.Task)
    def task(self, task):
        """
        Update the running statistics with a task result
        """

        t = systime.time()

        component = 'TASK'

        self.logger.update (t, component, 'host'                  , task.host)
        self.logger.update (t, component, 'tag'                   , task.tag)
        self.logger.update (t, component, 'result'                , task.result)
        self.logger.update (t, component, 'return_status'         , task.return_status)
        self.logger.update (t, component, 'total_bytes_transfered', task.total_bytes_transferred)

        # try/except: support different versions of cctools
        try:
            ### autobuild
            comp_time = task.cmd_execution_time
            comp_name = 'cmd_execution_time'
            self.logger.update (t, component, 'time_send_files'  ,
                                (task.send_input_finish - task.send_input_start     ) / 10.**6)
            self.logger.update (t, component, 'time_receive_files' ,
                                (task.receive_output_finish - task.receive_output_start ) / 10.**6)



        except AttributeError:
            ### "stable" version
            comp_time  = task.computation_time
            comp_name  = 'computation_time'

        ### convert all times to seconds from microseconds
        self.logger.update (t, component, comp_name            ,
                            comp_time                                                 / 10.**6)
        self.logger.update (t, component, 'total_transfer_time',
                            task.total_transfer_time                                  / 10.**6)

        self.logger.update (t, component, 'turnaround_time'    ,
                            (task.finish_time           - task.submit_time          ) / 10.**6)



    def save(self, wqstats, taskstats):
        """
        @param wqstats: path to save the wq stats to
        @param taskstats: path to save the task stats to
        """

        with open(wqstats, 'w') as fd:
            self._save_wq_stats(fd)

        with open(taskstats, 'w') as fd:
            self._save_task_stats(fd)


    def _save_attrs(self, fd, name, times, attrs):
        print 'Saving', name, 'data to', fd.name
        data = dict()
        data['time'] = times
        print '\t', 'time'
        for a in attrs:
            print '\t', a
            data[a] = getattr(self, a).values
        np.savez(fd, **data)

    def _save_task_stats(self, fd):
        attrs = 'computation_time total_bytes_transferred total_transfer_time task_life_time'.split()

        self._save_attrs(fd, 'task', self._task_times.get(), attrs)

    def _save_wq_stats(self, fd):
        attrs =  'workers_ready workers_busy tasks_running tasks_waiting tasks_complete'.split()
        attrs += 'total_tasks_dispatched total_tasks_complete total_workers_joined'.split()
        attrs += 'total_workers_removed total_bytes_sent total_bytes_received'.split()
        attrs += 'total_send_time total_receive_time'.split()

        self._save_attrs(fd, 'wq', self._wq_times.get(), attrs)



class Timings(object):

    def __init__(self):
        self.timer = Timer()
        self.times = ExtendableArray()
        self.stats = Statistics()

    @property
    def data(self):
        return self.times.get(), self.stats.values

    def start(self):
        self.timer.start()

    def stop(self):
        self.timer.stop()
        self.times.append(systime.time())
        self.stats.update(self.timer.elapsed())



class AWEStats(object):

    def __init__(self, logger=None):

        self.iteration = Timer()
        self.resample  = Timer()
        self.barrier   = Timer()

        self.logger    = logger or StatsLogger()

    @typecheck(str, Timings)
    def _timeit(self, state, timings, name):
        """
        *state* = {start|stop}
        """

        if state.lower() == 'start':
            timings.start()
        elif state.lower() == 'stop':
            timings.stop()
            t = systime.time()
            self.logger.update(t, 'AWE', name, timings.elapsed())
            timings.reset()
        else:
            raise ValueError, 'Unknown state %s: valid: {start|stop}' % state

    def time_iter(self, state):
        self._timeit(state, self.iteration, 'iteration time')

    def time_resample(self, state):
        self._timeit(state, self.resample, 'resample time')

    def time_barrier(self, state):
        self._timeit(state, self.barrier, 'barrier time')

    def save(self, path):
        print 'Saving to', path
        with open(path, 'w') as fd:
            data = dict()

            print '\t','iteration data'
            ts, vs = self.iteration.data
            data['iteration_time'] = ts
            data['iteration_values'] = vs

            print '\t','resample data'
            ts, vs = self.resample.data
            data['resample_time'] = ts
            data['resample_values'] = vs

            print '\t','barrier data'
            ts, vs = self.barrier.data
            data['barrier_time'] = ts
            data['barrier_values'] = vs

            np.savez(fd, **data)


class StatsLogger (object):

    def __init__(self, path='stats.log.bz2', buffersize=9):

        print 'StatsLogger opening', path
        self._fd = bz2.BZ2File(path, 'w', buffersize)
        self._path = path

    def __del__(self):
        self._fd.close()

    @property
    def path(self): return self._path

    @typecheck(float, str, str)
    def update(self, t, component, name, val):
        s = '%f %s %s %s\n' % (t, component, name, val)
        self._fd.write(s)

    def output(self, val):
        self._fd.write(str(val))
