
import awe


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

        list0          = np.arange(len(walkers)) ## walker indices
        weights        = walkers.weights         ## walker weights
        ntargetwalkers = self.targetwalkers      ## target number of walkers


        list1          = list()                  ## new list of walkers
        newweights     = list()                  ## weights of new walkers
        nwalkerlist1   = 0                       ## number of walkers in list 1

        for cell in set(walkers.cells):

            awe.log('Processing cell %s' % cell)

            ### initialize the list of weights for walkers in the current cell
            ixs = np.where(walkers.cells == cell)
            wi  = weights[ixs]
            ind = argsort(-wi)
            awe.log('\tsorted indexes %s' % ind)

            ### sort the walkers in descending order based on their weights
            awe.log('\tlist0: %s' % list0)
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
