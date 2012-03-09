
from util import typecheck
import aweclasses
import mdtools

import numpy as np
import itertools
import time


class IResampler(object):

    """
    Interface that all resampling methods should implement
    """

    def resample(self, walkers):
        """
        Parameters:
          *walkers* : an object of type awe.aweclasses.WalkerGroup

        Returns:
          a new awe.aweclasses.WalkerGroup instance
        """

        raise NotImplementedError

    @typecheck(aweclasses.WalkerGroup)
    def __call__(self, walkers):
        print time.asctime(), 'Resampling'
        ws2 = self.resample(walkers)
        assert type(ws2) is aweclasses.WalkerGroup
        return ws2



class Identity(IResampler):

    """
    the identity function
    """

    def resample(self, walkers):
        return walkers


class OneColor(IResampler):

    """
    A single/no color algorithm based on Eric Darve and Ernest Ryu's:
      "Computing reaction rates in bio-molecular systems using discrete macro-states"
    """

    def __init__(self, targetwalkers):
        self.targetwalkers = targetwalkers

    def resample(self, walkergroup):

        from numpy import floor, argsort, random, sum

        newwalkers = list()

        for cell in set(walkergroup.cells):
            # print 'Processing cell', cell

            ### initialize the list of weights for walkers in the current cell
            ixs        = np.where(walkergroup.cells == cell)
            oldwalkers = walkergroup.getslice(ixs)
            weights    = oldwalkers.weights.copy()
            # print '\tixs:', ixs
            # print '\tweights:', weights

            ### sort the walkers in descending order based on their weights,
            ##+ this ensures only walkers whose weight > targetWeight are split.
            mywalkers = list(argsort(-weights))
            # print '\tmywalkers:', mywalkers

            ### sanity check
            testmaxw = float('inf')
            for i in mywalkers:
                myw = weights[i]
                assert myw == oldwalkers[i].weight, 'Weights mismatch'
                assert myw <= testmaxw, 'Weights non-monotonically decreasing'
                testmaxw = myw

            ### setup cell weight and target weights
            W     = sum(weights)
            tw    = W / self.targetwalkers
            # print '\tW', W, 'tw', tw


            ### we assume that there is at least one walker in the cell
            x = mywalkers.pop()

            ### keep track of active walkers for splitting
            activewalkers = 0

            ### The algorithm terminates since the last walker removed
            ##+ from 'mywalkers' when W == tw. The max number of
            ##+ iterations is bounded by 'len(where(group.cell = cell) + targetwalkers'
            while True: # exit using break

                Wx = weights[x]
                currentWalker = oldwalkers[x]
                # print '\tweight of', x, 'is', Wx

                ### split
                if Wx > tw or len(mywalkers) == 0:

                    ### choose number of times to split
                    ##+ r = floor( Wx / tw )
                    ### min, max: work around round-off errors
                    r = max(1, int(floor( Wx/tw )) )
                    r = min(r, self.targetwalkers - activewalkers)
                    activewalkers += r
                    # print '\tactive walkers', activewalkers

                    ### split the current walker
                    # print '\tsplitting', x, r, 'times'
                    for _ in itertools.repeat(x, r):
                        w = aweclasses.Walker(start  = currentWalker.end,
                                              weight = tw,
                                              color  = currentWalker.color,
                                              cell   = cell)
                        newwalkers.append(w)


                    ### update the weights for the current walker and mark
                    ##+ for reconsideration
                    if activewalkers < self.targetwalkers and Wx - r * tw > 0.0:
                        mywalkers.append(x)
                        weights[x] = Wx - r * tw
                        # print '\tupdated weights of', x

                    ### continue the loop?
                    if len(mywalkers) > 0:
                        x = mywalkers.pop()
                    else: break

                ### merge
                else:
                    y = mywalkers.pop()
                    # print '\tmerging', x, y
                    Wy = weights[y]
                    Wxy = Wx + Wy
                    p = np.random.random()
                    if p < Wy / Wxy:
                        x = y
                    weights[x] = Wxy


        ### setup the WalkerGroup to return
        newgroup = aweclasses.WalkerGroup(count=len(newwalkers), topology=walkergroup.topology)
        for w in newwalkers:
            newgroup.add(w)

        return newgroup


class IPlotter(IResampler):

    def __init__(self, **kws):
        self.plotfile = kws.pop('plotfile', 'plot.png')

    def compute(self, walkergroup):
        raise NotImplementedError

    def plot(self):
        raise NotImplementedError

    def __call__(self, walkers):
        ws = IResampler.__call__(self, walkers)
        self.compute(ws)
        self.plot()
        return ws


class OneColor_SaveWeights(OneColor):

    def __init__(self, nwalkers, datfile='weights.dat'):
        OneColor.__init__(self, nwalkers)
        self.datfile   = datfile
        self.iteration = 0

    def saveweights(self, group, mode='a'):
        print 'Saving weights to', self.datfile

        ### all the walkers in a cell have the same weight, so we only
        ### need to save the (iteration, cell, weight) triples
        cells   = np.array(list(set(group.cells)))
        iters   = self.iteration * np.ones(len(cells))
        weights = -1 * np.ones(len(cells))
        for i, c in enumerate(cells):
            ixs        = np.where(group.cells == c)
            walkers    = group.getslice(ixs)
            w          = walkers.weights[0] ### assume at least one walker per cell
            weights[i] = w
        assert weights.min() >= 0
        assert weights.max() <= 1
        vals = np.vstack( (iters, cells, weights) )

        with open(self.datfile, mode) as fd:
            np.savetxt(fd, vals.T)

    def resample(self, walkergroup):
        if self.iteration == 0:
            self.saveweights(walkergroup, mode='w')

        newgroup         = OneColor.resample(self, walkergroup)
        self.iteration  += 1
        self.saveweights(newgroup)

        return newgroup
