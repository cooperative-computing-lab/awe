
import awe
import mdtools


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

    """
    Implementation based on Eric Darve and Ernest Ryu's paper:
    "Computing reaction rates in bio-molecular systems using discrete macro-states"
    """

    def __init__(self, targetwalkers):
        self.targetwalkers = targetwalkers

    def resample(self, walkers):

        awe.log('RESAMPLE weights: %s' % walkers.weights)

        from   numpy import floor, argsort, random, sum
        import numpy as np

        for cell in set(walkers.cells):

            list0          = np.arange(len(walkers)) ## walker indices
            weights        = walkers.weights         ## walker weights
            ntargetwalkers = self.targetwalkers      ## target number of walkers

            list1          = list()                  ## new list of walkers
            newweights     = list()                  ## weights of new walkers
            nwalkerlist1   = 0                       ## number of walkers in list 1


            awe.log('Processing cell %s' % cell)

            ### initialize the list of weights for walkers in the current cell
            ixs = np.where(walkers.cells == cell)
            wi  = weights[ixs]
            ind = argsort(-wi)

            ### sort the walkers in descending order based on their weights
            list0 = list(list0[ind])

            W     = sum(wi)
            tw    = W / ntargetwalkers

            ### we assume that there is at least on walker in the cell
            x     = list0.pop()

            while True: ## exit using a break

                Wx = weights[x]
                if (Wx > tw or len(list0) == 0):

                    ## min, max required because of round-off errors
                    print 'Wx', Wx, 'tw', tw
                    r = max(1, int( floor(Wx/tw) ))
                    r = min(r, ntargetwalkers - nwalkerlist1)

                    ### update the number of walkers in list1
                    nwalkerlist1 += r

                    ### insert r copies of walkers in list1
                    for item in xrange(x, r):
                        list1.append(item)
                        newweights.append(tw)

                    if nwalkerlist1 < ntargetwalkers and Wx - r*tw > 0.0:
                        list0.append(x)
                        weights[x] = Wx - r*tw

                    if len(list0) > 0:
                        x = list0.pop()
                    else:
                        break

                else:

                    y   = list0.pop()
                    Wy  = weights[y]
                    Wxy = Wx + Wy

                    ## randomly select a walker
                    p   = random.random()
                    if p < Wy / Wxy:
                        x = y

                    weights[x] = Wxy

        walkers.weights = weights
        return walkers


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
            plt.scatter(xs, ws)

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
