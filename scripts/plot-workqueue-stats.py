#!/usr/bin/env python

import pylab as plt
import numpy as np

import collections
import os, sys, argparse

WQStats = collections.namedtuple(
    'WQStats',
    ['timestamp',
     'start_time',
     'workers_init',
     'workers_ready',
     'workers_busy',
     'workers_cancelling',
     'tasks_waiting',
     'tasks_running',
     'tasks_complete',
     'total_tasks_dispatched',
     'total_tasks_complete',
     'total_workers_joined',
     'total_workers_connected',
     'total_workers_removed',
     'total_bytes_sent',
     'total_bytes_received',
     'total_send_time',
     'total_receive_time',
     'efficiency',
     'idle_percentage',
     'capacity',
     'avg_capacity',
     'port',
     'priority',
     ])


PLOT_ATTRS  = 'workers_init workers_ready workers_busy workers_cancelling'.split()
PLOT_ATTRS += 'tasks_waiting tasks_running tasks_complete'.split()
PLOT_ATTRS += 'total_tasks_dispatched total_tasks_complete'.split()
PLOT_ATTRS += 'total_workers_joined total_workers_connected total_workers_removed'.split()
PLOT_ATTRS += 'total_bytes_sent total_bytes_received'.split()
PLOT_ATTRS += 'total_send_time total_receive_time'.split()
PLOT_ATTRS += 'efficiency idle_percentage capacity avg_capacity'.split()



def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-i','--logpath', help='Path to WQ logfile')
    p.add_argument('-o', '--plotpath', help='Where to store the figure')

    p.add_argument('-a', '--plot-all', action='store_true', help='Plot all attributes')

    for a in PLOT_ATTRS:
        p.add_argument('--' + a, action='store_true')

    return p.parse_args()

def load_wqstats(path):
    records = np.loadtxt(path, unpack=True)
    stats   = WQStats(*list(records))
    return stats

def adjust(wqs):
    # convert to seconds
    s = wqs._replace(timestamp = (wqs.timestamp - wqs.timestamp[0]) / 10.**6)

    return s

def plot(opts, wqs):

    for a in PLOT_ATTRS:
        if opts.plot_all or getattr(opts, a):
            print 'Plotting', a
            plt.plot(wqs.timestamp, getattr(wqs, a), label=a)

    plt.legend(loc='upper left')

    plt.xlabel('Time (s)')

    plt.savefig(opts.plotpath)

if __name__ == '__main__':
    opts = getopts()
    wqstats = load_wqstats(opts.logpath)
    wqstats = adjust(wqstats)
    plot(opts, wqstats)
