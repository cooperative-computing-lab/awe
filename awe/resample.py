
import awe
import mdtools

import numpy as np
import itertools


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

    @awe.trace()
    def __call__(self, walkers):
        assert type(walkers) is awe.aweclasses.WalkerGroup
        ws2 = self.resample(walkers)
        assert type(ws2) is awe.aweclasses.WalkerGroup
        return ws2



class Identity(IResampler):

    """
    the identity function
    """

    @awe.trace()
    def resample(self, walkers):
        return walkers


class OneColor(IResampler):

    def __init__(self, targetwalkers):
        self.targetwalkers = targetwalkers

        ## required because of roundoff errors
        self.eps           = np.finfo(np.double).eps

    def resample(self, walkergroup):

        print 'Resampling'

        from numpy import floor, argsort, random, sum

        newwalkers = list()

        for cell in set(walkergroup.cells):
            print 'Processing cell', cell

            ### initialize the list of weights for walkers in the current cell
            ixs        = np.where(walkergroup.cells == cell)
            oldwalkers = walkergroup.getslice(ixs)
            weights    = oldwalkers.weights.copy()
            print '\tixs:', ixs
            print '\tweights:', weights

            ### sort the walkers in descending order based on their weights,
            ##+ this ensures only walkers whose weight > targetWeight are split.
            mywalkers = list(argsort(-weights))
            print '\tmywalkers:', mywalkers

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
            print '\tW', W, 'tw', tw


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
                print '\tweight of', x, 'is', Wx

                ### split
                if Wx > tw or len(mywalkers) == 0:
                    ### work around round-off errors
                    r = max(1, int(floor( Wx/tw )) )
                    r = min(r, self.targetwalkers - activewalkers)
                    activewalkers += r
                    print '\tactive walkers', activewalkers

                    ### split the current walker
                    print '\tsplitting', x, r, 'times'
                    for _ in itertools.repeat(x, r):
                        w = awe.aweclasses.Walker(coords = currentWalker.coords,
                                                  weight = tw,
                                                  color  = currentWalker.color,
                                                  cell   = cell)
                        newwalkers.append(w)


                    ### update the weights for the current walker and mark
                    ##+ for reconsideration
                    if activewalkers < self.targetwalkers and Wx - r * tw > 0.0:
                        mywalkers.append(x)
                        weights[x] = Wx - r * tw
                        print '\tupdated weights of', x

                    ### continue the loop?
                    if len(mywalkers) > 0:
                        x = mywalkers.pop()
                    else: break

                ### merge
                else:
                    y = mywalkers.pop()
                    print '\tmerging', x, y
                    Wy = weights[y]
                    Wxy = Wx + Wy
                    p = np.random.random()
                    if p < Wy / Wxy:
                        x = y
                    weights[x] = Wxy


        newgroup = awe.aweclasses.WalkerGroup(count=len(newwalkers), topology=walkergroup.topology)
        for w in newwalkers:
            newgroup.add(w)

        ### normalize weights
        newgroup.weights /= np.sum(newgroup.weights)
        return newgroup


class SimpleWeightsPlotter(Simple):
    def __init__(self, nwalkers, plotfile='weights.png'):
        Simple.__init__(self, nwalkers)
        self._plotfile = plotfile
        self._weights  = list()

    @property
    def plotfile(self):
        return self._plotfile

    def __getitem__(self, i):
        return self._weights[i]

    def __len__(self):
        return len(self._weights)


    def resample(self, walkers):
        self._weights.append(walkers.weights.copy())
        self.plot()
        return Simple.resample(self, walkers)

    def plot(self):
        print 'Plotting'
        import matplotlib.pylab as plt
        import numpy as np

        for i in xrange(len(self)):
            ws = self[i]
            xs = i * np.ones(len(ws), dtype=int)
            for x, w in zip(xs, ws):
                plt.scatter(x, w, alpha=0.75, color=plt.cm.jet(1.*i/len(xs)))

        plt.ylabel('Weights')
        plt.xlabel('Iteration')

        plt.savefig(self.plotfile)
        print 'Saved image to', self.plotfile


class SimpleRMSDPlotter(Simple):
    def __init__(self, nwalkers, ref, sel='name CA', plotfile='rmsds.png'):

        Simple.__init__(self, nwalkers)
        self.plotfile = plotfile
        self.rmsds    = list()
        self.sel      = sel
        self.ref      = mdtools.pdb.load(ref)
        self.ref.setunits('nanometers')

    def resample(self, walkers):
        coords = walkers.positions
        ref    = self.ref.atomcoords

        print 'Computing RMSD'
        rmsd = mdtools.analysis.Analysis.rmsd(ref, *coords)
        self.rmsds.append(rmsd)
        print self.rmsds
        self.plot()

        print 'Resampling'
        return Simple.resample(self, walkers)


    def plot(self):
        print 'Plotting'
        import matplotlib.pylab as plt
        import numpy as np

        plt.jet()

        for i in xrange(len(self.rmsds)):
            rs = self.rmsds[i]
            xs = np.arange(len(rs), dtype=int)
            plt.scatter(xs, rs, alpha=0.5, color=plt.cm.jet(1.*i/len(self.rmsds)))
        plt.xlabel('Walker')
        plt.ylabel('RMSD')

        print 'Saving plot to', self.plotfile
        plt.savefig(self.plotfile)
