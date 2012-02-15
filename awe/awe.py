

import awe

import mdtools

import numpy as np




class AWE(object):

    def __init__(self, wqconfig=None, cells=None, walkers=None, iterations=-1, resample=None):

        assert type(wqconfig) is awe.workqueue.Config
        # TODO: assert type(cells) is
        assert type(walkers) is WalkerGroup
        assert type(iterations) is int
        # TODO: assert type(resample) is

        self.wq         = awe.workqueue.WorkQueue(wqconfig)
        self.cells      = cells
        self.walkers    = walkers
        self.iterations = iterations
        self.resample   = resample

        self.stats      = awe.stats.AWEStats()


    def _submit(self):

        for i in xrange(len(self.walkers)):
            task = self.wq.new_task(self.walkers, i)

            self.wq.submit(task)


    def _recv(self):

        self.stats.time_barrier('start')
        while not self.wq.empty:
            awe.log('Waiting for result')
            walker     = self.wq.recv()
            walkers[k] = walker
        self.stats.time_barrier('stop')


    def _resample(self):

        self.stats.time_resample('start')
        self.walkers = self.resample(self.walkers)
        self.stats.time_resample('stop')
            

    def run(self):

        for iteration in xrange(self.iterations):

            self.stats.time_itr('start')

            self._submit()
            self._recv()     ## barrier
            self._resample()

            self.stats.time_itr('stop')
