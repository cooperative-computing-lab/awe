
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


class Simple(IResampler):

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
            ixs       = np.where(walkergroup.cells == cell)
            ixs       = ixs[0]
            weights   = walkergroup.weights[ixs]
            print '\tixs:', ixs

            ### sort the walkers in descending order based on their weights
            ##+ this ensures only walkers whose weight > targetWeight are split
            ixs = ixs[argsort(-weights)]
            mywalkers = list(ixs)
            weights   = walkergroup.weights.copy()

            print '\tmywalkers:', mywalkers
            print '\tweights:', weights[mywalkers]
            

            ### setup total weight and target weights
            W     = sum(weights)
            tw    = W / self.targetwalkers
            print '\tW', W, 'tw', tw

            ### TODO
            activewalk = 0

            if len(mywalkers) > 0:

                ### x,y are walker indices
                x = mywalkers.pop()
                print'\twalker x:', x
                currentWalker = walkergroup[x]

                while True:

                    Wx = weights[x]

                    ### split
                    if (Wx + self.eps > tw):
                        print 'a'
                        ### determine number copies of current walker needed
                        r = int(np.floor( (Wx+self.eps) / tw ))

                        ### split: insert r copies of walkers in list1
                        for item in itertools.repeat(x, r):
                            print 'b'
                            w = awe.aweclasses.Walker(coords = currentWalker.coords,
                                                      weight = tw,
                                                      color  = currentWalker.color,
                                                      cell   = cell)
                            newwalkers.append(w)

                        ### ???
                        activewalk += r
                        if activewalk < self.targetwalkers and Wx-r*tw+self.eps > 0.0:
                            print 'c'
                            w = awe.aweclasses.Walker(coords = currentWalker.coords,
                                                      weight = Wx - r * tw,
                                                      color  = currentWalker.color,
                                                      cell   = cell)
                            mywalkers.append(x)
                            newwalkers.append(w)

                        if len(mywalkers) > 0:
                            print 'd'
                            x = mywalkers.pop()
                        else: break

                    ### merge
                    else:
                        print 'e'
                        if len(mywalkers) > 0:
                            print 'f'
                            y = mywalkers.pop()
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
            plt.scatter(xs, ws, alpha=0.75, color=plt.cm.jet(1.*i/len(self)))

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
