"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


from util import typecheck, returns
import aweclasses

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

    @typecheck(aweclasses.System)
    @returns(aweclasses.System)
    def __call__(self, s1):
        print time.asctime(), 'Resampling'
        return self.resample(s1)


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

    def resample(self, system):

        from numpy import floor, argsort, random, sum

        newsystem = system.clone()

        for cell in system.cells:

            ### initialize the list of weights for walkers in the current cell
            localsystem = system.filter_by_cell(cell)
            weights    = localsystem.weights
            walkers    = localsystem.walkers

            print time.asctime(), 'Resampling cell', cell, len(walkers) ,'walkers'

            if not len(walkers) > 0: continue

            ### sort the walkers in descending order based on their weights,
            ##+ this ensures only walkers whose weight > targetWeight are split.
            mywalkers = list(argsort(-weights))
            print '\tmywalkers:', mywalkers

            ### sanity check
            testmaxw = float('inf')
            for i in mywalkers:
                myw = weights[i]
                assert myw == weights[i], 'Weights mismatch'
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
            ##+ iterations is bounded by 'len(where(system.cell == cell) + targetwalkers'
            while True: # exit using break

                Wx = weights[x]
                currentWalker = walkers[x]
                print '\tweight of', x, 'is', Wx

                ### split
                if Wx >= tw or len(mywalkers) == 0:

                    ### choose number of times to split
                    ##+ r = floor( Wx / tw )
                    ### min, max: work around round-off errors
                    r = max(1, int(floor( Wx/tw )) )
                    r = min(r, self.targetwalkers - activewalkers)
                    activewalkers += r
                    print '\tactive walkers', activewalkers

                    ### split the current walker
                    print '\tsplitting', x, r, 'times'
                    for _ in itertools.repeat(x, r):
                        w = currentWalker.restart(weight=tw)
                        newsystem.add_walker(w)


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

        return newsystem

class MultiColor(OneColor):

    def __init__(self, nwalkers, partition):
        OneColor.__init__(self, nwalkers)
        self.partition   = partition
        ncolors          = partition.ncolors
        self.transitions = np.zeros((ncolors, ncolors))

    def resample(self, system):

        ### update colors
        for w in system.walkers:
            cell     = system.cell(w.assignment)

            oldcolor = w.color
            newcolor = self.partition.color(cell)
            if newcolor is None: newcolor = oldcolor

            if not oldcolor == newcolor:
                print 'Updating color:', w, oldcolor, '->', newcolor
                w.color = newcolor

            self.transitions[oldcolor, newcolor] += w.weight

        ### resample individual colors using OneColor algorithm
        newsystem = aweclasses.System(topology=system.topology)
        for color in system.colors:
            thiscolor  = system.filter_by_color(color)
            print time.asctime(), 'Resampling color', color, len(thiscolor.walkers), 'walkers'
            resampled  = OneColor.resample(self, thiscolor)
            newsystem += resampled


        return newsystem

    def save_transitions(self, path):
        print time.asctime(), 'Saving transition matrix to', repr(path)
        np.savetxt(path, self.transitions)


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

class ISaver(IResampler):

    """
    Save after each resampling procedure
    """

    @typecheck(IResampler, datfile=str)
    def __init__(self, resampler, datfile='save.dat'):
        self.resampler = resampler
        self.datfile   = datfile
        self.iteration = 0

    @typecheck(aweclasses.System, mode=str)
    def save(self, system, mode='a'):
        self._save(system, mode=mode)

    @returns(str)
    def heading(self):
        return self._heading()

    def _save(self, system, mode='a'):
        raise NotImplementedError

    def _heading(self):
        return ''

    @typecheck(aweclasses.System)
    @returns(aweclasses.System)
    def resample(self, system):
        if self.iteration == 0:
            with open(self.datfile, 'w') as fd:
                fd.write(self.heading())

        newsystem = self.resampler.resample(system)
        self.iteration += 1
        self.save(newsystem, mode='a')

        return newsystem

class SaveWeights(ISaver):

    def __init__(self, resampler, datfile='weights.dat'):
        ISaver.__init__(self, resampler, datfile=datfile)

    def _heading(self):
        return \
            '# Each line represents a walker at:\n' + \
            '# walkerid iteration cell weight color\n'

    def _save(self, system, mode='a'):
        print time.asctime(), 'Saving weights to', self.datfile

        ### all the walkers in a cell have the same weight, so we only
        ### need to save the walkerid, iteration, cell, weight, and color for each walker

        with open(self.datfile, mode) as fd:
            for i, w in enumerate(system.walkers):
                s = '%(wid)d\t%(iteration)d\t%(cell)d\t%(weight)f\t%(color)s\n' % {
                    'wid'       : w.id           ,
                    'iteration' : self.iteration ,
                    'cell'      : w.assignment   ,
                    'weight'    : w.weight       ,
                    'color'     : w.color } # if w.color is not None else 'nan'       }
                fd.write(s)
