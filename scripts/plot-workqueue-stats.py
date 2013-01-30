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

    p.add_argument('-H', '--history', type=float, help='Number of hours to plot [default: all]')

    p.add_argument('-a', '--plot-all', action='store_true', help='Plot all attributes')

    for a in PLOT_ATTRS:
        p.add_argument('--' + a, action='store_true')

    return p.parse_args()

def load_wqstats(path):
    print 'Loading data from', path
    records = np.loadtxt(path, unpack=True)
    stats   = WQStats(*list(records))
    return stats

def adjust(wqs, history=None):
    s = wqs

    print 'Attempting to adjust history window'
    s = adjust_history(s, history=history)

    print 'Converting timestamp to hours'
    s = s._replace(timestamp = s.timestamp / (3600. * 10.**6))

    print 'Normalizing to start time'
    s = s._replace(timestamp = s.timestamp - s.timestamp[0])

    return s

def adjust_history(wqs, history=None):
    if history is None:
        return wqs
    else:
        now = wqs.timestamp[-1]
        hus = history * 3600 * 10**6  # convert from hours
        ago = now - hus
        ixs = np.where(wqs.timestamp >= ago)
        print 'Filtered',  len(wqs.timestamp) - len(wqs.timestamp[ixs]), 'values'

        d = wqs._asdict()
        for a, v in d.items():
            d[a] = v[ixs]

        return WQStats(**d)


def plot(opts, wqs):

    for a in PLOT_ATTRS:

        if opts.plot_all or getattr(opts, a):
            print 'Plotting', a
            v = getattr(wqs, a)
            plt.plot(wqs.timestamp, v, label=a)

    plt.legend(loc='upper left')

    plt.xlabel('Time (hours)')

    plt.savefig(opts.plotpath)

if __name__ == '__main__':
    opts = getopts()
    wqstats = load_wqstats(opts.logpath)
    wqstats = adjust(wqstats, history=opts.history)

    plot(opts, wqstats)
