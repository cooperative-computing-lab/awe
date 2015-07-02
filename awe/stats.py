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
import gzip
import os


class Timer(object):
    
    """
    awe.stats.Timer

    A stopwatch-esque functionality that uses system time.

    Fields:
        t0 - Start time
        t1 - End time

    Methods:
        reset       - erase the contents of t0 and t1
        start       - set t0 to the current time
        stop        - set t1 to the current time
        isrunning   - determine if the stopwatch is currently recording
        elapsed     - output the time elapsed between t0 and t1 or the current time
    """
    
    def __init__(self):
        """
        awe.stats.Timer.__init__

        Alias for Timer.reset called upon instance creation.

        Parameters:
            self - the Timer instance accessed

        Returns:
            None
        """
        
        self.reset()


    def reset(self):
        """
        awe.stats.Timer.reset

        Clears any stored times and sets t0 and t1 to native states.

        Parameters:
            self - the Timer instance accessed

        Returns:
            None
        """
        
        self.t0 = 0
        self.t1 = float('inf')

    def start(self):
        """
        awe.stats.Timer.start

        Sets t0 to the current time, "starting" the stopwatch.

        Parameters:
            self - the Timer instance accessed

        Returns:
            None
        """
        
        self.t0 = systime.time()

    def stop(self):
        """
        awe.stats.Timer.stop

        Sets t1 to the current time, "stopping" the stopwatch.

        Parameters:
            self - the Timer instance accessed

        Returns:
            None
        """
        
        if self.t0 > 0:
            self.t1 = systime.time()

    def isrunning(self):
        """
        awe.stats.Timer.isrunning

        Determines whether the stopwatch is active or not, in the sense that
        a start time has been recorded and a stop time has not.

        Parameters:
            self - the Timer instance accessed

        Returns:
            A Boolean value indicating whether or not the stopwatch is active
        """
        
        return self.t0 > 0 and self.t1 == float('inf')

    def elapsed(self, current=True, units='s'):
        """
        awe.stats.Timer.elapsed

        Returns the elapsed time since t0 started if the timer is running.

        Parameters:
            self    - the Timer instance being accessed
            current - whether or not to use the current time as a reference
            units   - the time unit to which the elapsed time should be
                      converted ['s' = seconds (default), 'm' = minutes,
                      'h' = hours, 'd' = days]

        Returns:
            A number representing the elapsed time in the specified units. If
            the Timer instance has never been started, returns 0.
        """
        
        if self.t0 > 0:
            # System time is measured in seconds, so convert from seconds
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

            # diff is the elapsed time in seconds
            diff       = t1 - self.t0

            # Convert to the specified units
            return multiplier * diff
        
        # Only reached if the Timer instance was not started 
        return 0


### use a single global timer for the system
_TIMER = Timer()

### simulate built-in *time* module
class time:
    
    """
    awe.stats.time

    Simulates the system time module using a singleton of awe.stats.Timer as a
    starting point instead of the UNIX epoch.

    Fields:
        None

    Methods:
        start   - (static) start the Timer singleton
        time    - (static) return the elapsed time from the Timer singleton
        timer   - (static) return a reference to the Timer singleton
    """
    
    @staticmethod
    def start():
        """
        awe.stats.time.start

        Starts the Timer singleton.

        Parameters:
            None

        Returns:
            None
        """
        
        global _TIMER
        _TIMER.start()

    @staticmethod
    def time():
        """
        awe.stats.time.time

        Determines the elapsed time of the Timer singleton.

        Parameters:
            None

        Returns:
            A number representing the elapsed time in seconds since the Timer
            singleton was started
        """
        
        global _TIMER
        t = _TIMER.elapsed(units='s')
        return t


    @staticmethod
    def timer():
        """
        awe.stats.time.timer

        Allows direct access to the Timer singleton.

        Parameters:
            None

        Returns:
            A reference to the Timer singleton object.
        """
        
        global _TIMER
        return _TIMER


class ExtendableArray(object):
    
    """
    awe.stats.ExtendableArray

    Contains and manages an oversized numpy.array object for efficient append
    operations.

    Fields:
        _type   - a numpy array type (e.g. numpy.zeros, numpy.ones, etc.) by
                  which to create the numpy.array
        _size0  - the total capacity of the numpy.array
        _count  - a counter keeping track of the number of entries in the array
        _factor - the multiplicative factor by which to increase the array
                  capacity
        _vals   - the numpy.array object

    Methods:
        _initialize     - create a new numpy.array
        get             - return the contents of the array up to the current
                          count
        _realloc        - allocate a new numpy.array with increased capacity
        append          - add elements to the end of the numpy.array
        __getitem__     - retrieve an element from the numpy.array
        __setitem__     - change the value of an existing element
        __len__         - return the length (count, not capacity) of the array
        __contains__    - determine if the array contains a specified element
        __iter__        - return an iterator over the set values of the array
        __str__         - return a string representation of the array
        __repr__        - return a string representation of the array
    """

    def __init__(self, typ=np.zeros, size=500, factor=2):
        """
        ExtendableArray.__init__

        Set instance fields at initialization.

        Parameters:
            typ     - a numpy array type (e.g. numpy.zeros, numpy.ones, etc.) by
                      which to create the numpy.array
            size    - the initial capacity of the numpy.array
            factor  - the multiplicative factor by which to increase the array
                      capacity
        """
        self._type    = typ
        self._size0   = size
        self._count   = 0       # Says no elements have been appended
        self._factor  = factor
        self._vals    = self._initialize()

    def _initialize(self):
        """
        awe.stats.ExtendableArray._initialize

        Create a new array of the ExtendableArray instance numpy.array type
        with capacity specified by the ExtendableArray instance.

        Parameters:
            self - the ExtendableArray instance accessed

        Returns:
            A new numpy.array of type self._type and length self._size0
        """
        
        return self._type(self._size0)

    def get(self):
        """
        awe.stats.ExtendableArray.get

        Return only the modified entries of the numpy.array managed by the 
        ExtendableArray instance.

        Parameters:
            self - the ExtendableArray instance accessed

        Returns:
            A slice of the numpy.array size self._count consisting only of
            elements appended to this instance of ExtendableArray.
        """
        
        return self._vals[:self._count]

    def _realloc(self):
        """
        awe.stats.ExtendableArray._realloc

        Resize the numpy.array managed by the instance of ExtendableArray if
        the capacity has been reached.
        
        Parameters:
            self - the ExtendableArray instance accessed

        Returns:
            None
        """
        
        if self._count == len(self._vals):

            x     = self._count
            alloc = x * self._factor
            self._size0 = alloc
            vals2 = self._initialize()
            vals2[:self._count] = self._vals[:self._count]
            self._vals = vals2


    def append(self, *values):
        """
        awe.stats.ExtendableArray.append

        Add a new value to the numpy.array managed by the ExtendableArray
        instance at the index of the current counter. Resize the array if
        capacity has been reached.
        
        Parameters:
            self    - the ExtendableArray instance accessed
            *values - the elements to append

        Returns:
            None
        """

        self._realloc() # Dangerous: What if len(values) >= _factor * _size0 ?
        i = self._count
        j = i + len(values)
        self._vals[i:j] = np.array(values) # Appends len(values) values
        self._count += len(values)

    def __getitem__(self, i):
        """
        awe.stats.ExtendableArray.__getitem__

        Retrieve the element at the specified index or exit if k is too large.
        The system exit seems like overkill. Maybe an Exception instead?
        
        Parameters:
            self    - the ExtendableArray instance accessed
            i       - the index from which to retrieve the element

        Returns:
            The element at index i if i is in range
        """

        if i < 0: k = self._count + i # Can look up in reverse. Nice!
        else:     k = i

        assert k <= self._count # Seems like overkill
        return self._vals[k]

    def __setitem__(self, i, v):
        """
        awe.stats.ExtendableArray.__setitem__

        Update the entry at the specified index with the supplied value. Once
        again, there's a strange assert here that seems like overkill.
        
        Parameters:
            self    - the ExtendableArray instance accessed
            i       - the index to access
            v       - the the new value

        Returns:
            None
        """

        if i < 0: k = self._count + i # Once again, reverse lookup
        else:     k = i

        assert k <= self._count # Seems like overkill
        self._vals[k] = v

    def __len__(self):
        """
        awe.stats.ExtendableArray.__len__

        Return the size of the active array.
        
        Parameters:
            self    - the ExtendableArray instance accessed

        Returns:
            The value of self._count, reflecting the size of the active array
        """

        return self._count

    def __contains__(self, v):
        """
        awe.stats.ExtendableArray.__contains__

        Determine whether or not a specified value is an entry in the array.
        
        Parameters:
            self    - the ExtendableArray instance accessed
            v       - the element to search for

        Returns:
            A Boolean value reflecting whether or not the supplied value is
            in the active array.
        """

        return v in self._vals[:self._count]

    def __iter__(self):
        """
        awe.stats.ExtendableArray.__iter__

        Return an iterator over the active array.
        
        Parameters:
            self    - the ExtendableArray instance accessed

        Returns:
            An iterator spanning the active array.
        """

        return iter(self._vals[:self._count])

    def __str__(self):
        """
        awe.stats.ExtendableArray.__str__

        Return a string visualization of the array contents.
        
        Parameters:
            self    - the ExtendableArray instance accessed

        Returns:
            A string representing the array contents
        """

        return str(self._vals[:self._count])

    def __repr__(self):
        """
        awe.stats.ExtendableArray.append

        Return a string representing the objects in the array.
        
        Parameters:
            self    - the ExtendableArray instance accessed

        Returns:
            A string representing the objects in the array
        """

        return repr(self._vals[:self._count])
            

class Statistics(object):

    """
    awe.stats.Statistics

    Keeps track of running statistics as workers work using an ExtendableArray.

    Fields:
        _num    - The number of of values in the ExtendableArray
        _mean   - The mean of values in the ExtendableArray
        _m2     - The second statistical moment
        _values - The ExtendableArray containing values on which to run stats

    Properties:
        num     - The number of elements in the ExtendableArray
        mean    - The mean of values in the ExtendableArray
        var     - The variance of values in the ExtendableArray
        std     - The standard deviation of values in the ExtendableArray

    Methods:
        values  - Getter for the ExtendableArray
        update  - Append values to the ExtendableArray and update the stats 

    """

    def __init__(self):

        self._num    = 0
        self._mean   = 0.
        self._m2     = 0.

        self._values = ExtendableArray() # keep track of the actual values for plots

    # Properties for useful statistics. The first two were probably unneccessary.
    num  = property(lambda self: self._num)
    mean = property(lambda self: self._mean)
    var  = property(lambda self: self._m2 / float(self._num))
    std  = property(lambda self: math.sqrt(self.var))

    @property
    def values(self):
        """
        awe.stats.Statistics.values

        Access the values added to the ExtendableArray.

        Parameters:
            None

        Returns:
            numpy.array slice containing all values added to the statistics class
        """

        return self._values.get()

    def update(self, *values):
        """
        awe.stats.Statistics.update

        Add a value to the ExtendableArray and update all statistics.
        
        Parameters:
            values - a list of values to add

        Returns:
            None
        """

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
    awe.stats.WQStats

    Keeps track of statistics about running time and the amount of data
    transferred.

    Fields:
        logger                  - the logging utility to use
        _task_times             - the list of running times for returned tasks
        _wq_times               - the list of running times for parts of WQ
        computation_time        - keeps track of statistics surrounding the
                                  computational time for a task
        total_bytes_transferred - keeps track of statistics surrounding the
                                  amount of data transferred
        total_transfer_time     - keeps track of statistics surrounding the
                                  amount of time spent transferring data
        task_life_time          - keeps track of statistics surrounding the
                                  amount of time tasks take to execute
    """

    def __init__(self, logger=None):

        """
        awe.stats.WQStats.__init__

        Create an instance of WQStats with a specified logging class.

        Parameters:
            logger - the logging class to use

        Returns:
            None
        """

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
        awe.stats.WQStats.task

        Add information about a returned task to the logger.

        Parameters:
            task - a completed Work Queue task

        Returns:
            None
        """

        t = systime.time()

        component = 'TASK'

        self.logger.update(t, component, 'host'                  , task.host)
        self.logger.update(t, component, 'tag'                   , task.tag)
        self.logger.update(t, component, 'result'                , task.result)
        self.logger.update(t, component, 'return_status'         , task.return_status)
        self.logger.update(t, component, 'total_bytes_transfered', task.total_bytes_transferred)

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

    def close(self):
        self.logger.close()

    def open(self):
        self.logger.open()


class StatsLogger (object):

    def __init__(self, path='debug/stats.log.gz', buffersize=9):

        prefix = os.path.dirname(os.path.abspath(os.path.expanduser(path)))
        if not os.path.exists(prefix):
            os.makedirs(prefix)

        self._fd = None
        self._path = path
        self._buffersize = buffersize

        self.open()

    @property
    def path(self): return self._path

    @typecheck(float, str, str)
    def update(self, t, component, name, val):
        s = '%f %s %s %s\n' % (t, component, name, val)
        self._fd.write(s)

    def output(self, val):
        self._fd.write(str(val))

    def close(self):
        if self._fd is not None:
            self._fd.close()
            self._fd = None

    def open(self):
        if self._fd is None:
            self._fd = gzip.GzipFile(self._path, 'ab')
