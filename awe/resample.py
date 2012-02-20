
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


class OneColor_PlotCellRMSD(OneColor,IPlotter):

    def __init__(self, nwalkers, **kws):
        OneColor.__init__(self, nwalkers)
        IPlotter.__init__(self, **kws)

        ref = kws.pop('ref')
        ref = mdtools.pdb.load(ref)
        ref.setunits('nanometers')

        self.ref   = ref.atomcoords
        self.rmsds = list()

    def compute(self, walkers):

        print 'Computing RMSDS'

        rmsds = list()
        for cell in set(walkers.cells):

            print '\tcell', cell

            ixs = np.where(walkers.cells == cell)
            coords = walkers.positions[ixs]
            rmsd = mdtools.analysis.Analysis.rmsd(self.ref, *coords)
            rmsd = np.mean(rmsd)
            rmsds.append((cell,rmsd))

        self.rmsds.append(rmsds)

    def plot(self):
        import matplotlib.pylab as plt

        print 'Plotting'

        for iteration, group in enumerate(self.rmsds):

            for cell, rmsd in group:
                print '\t', cell, rmsd
                plt.scatter(cell, rmsd, alpha=0.5, color=plt.cm.jet(1.*iteration/len(self.rmsds)))

        plt.xlabel('Cell')
        plt.ylabel('Average RMSD')
        plt.savefig(self.plotfile)

        print 'Saved figure to', self.plotfile

