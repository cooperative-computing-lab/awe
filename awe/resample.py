# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


from util import typecheck, returns, makedirs_parent
import aweclasses

import numpy as np
import itertools
import time, os
import textwrap

OUTPUT_DIR = os.getcwd()


class IResampler(object):

    """
    awe.resample.IResampler

    Interface that all resampling methods should implement.

    Fields:
        None

    Methods:
        resample - divide and kill off walkers until weights are equal to
                   the group's mean weight
    """

    def resample(self, walkers):
        """
        awe.IResampler.resample

        Placeholder for implementation in subclasses.

        Parameters:
          walkers - an awe.aweclasses.System instance

        Returns:
          A new awe.aweclasses.System instance
        """

        raise NotImplementedError

    @typecheck(aweclasses.System)
    @returns(aweclasses.System)
    def __call__(self, s1):
        print time.asctime(), 'Resampling'
        return self.resample(s1)


class Identity(IResampler):

    """
    awe.resample.Identity

    A resampler that does not alter the System it receives.

    Fields:
        None

    Methods:
        resample - returns the identity of the System
    """

    def resample(self, walkers):

        """
        awe.resample.Identity.resample

        Perform the identity operation on the System.

        Parameters:
            walkers - an awe.aweclasses.System instance

        Returns:
            The supplied awe.aweclasses.System instance
        """

        return walkers


class OneColor(IResampler):

    """
    awe.resample.OneColor

    A single/no color algorithm based on Eric Darve and Ernest Ryu's:
      "Computing reaction rates in bio-molecular systems using discrete
       macro-states"
    This algorithm does not try to determine which state the walker is in;
    instead, one state is assumed and all walkers are set to the same target
    weight.

    Fields:
        targetwalkers - a target weight that the new generation of walkers
                        should be set to
        histfile      - the file to which output weights are recorded

    Methods:
        resample - create a new generation of walkers by merging an splitting
                   weights from a processed generation
    """

    def __init__(self, targetwalkers):

        """
        awe.resample.OneColor.__init__

        Initialize a new instance of OneColor.

        Parameters:
            targetwalkers - a list of target weights to converge to (numeric)

        Returns:
            None
        """

        self.targetwalkers = targetwalkers
#####
        self.histfile = os.path.join(OUTPUT_DIR, 'walker-history.csv')
        makedirs_parent(self.histfile)
        with open(self.histfile, 'a') as fd:
            fd.write('%origID, parentID, currentID \n')

    def resample(self, system):

        """
        awe.resample.OneColor.resample

        Adjust weights of a processed group of walkers to meet the target group
        assuming there exists only one metastable state (macro-state).

        Parameters:
            system - an awe.aweclasses.System instance

        Returns:
            A new awe.aweclasses.System instance with adjusted weights
        """

        histfile_fd = open(self.histfile, 'a')

        from numpy import floor, argsort, random, sum

        newsystem = system.clone()

        for cell in system.cells:

            ### initialize the list of weights for walkers in the current cell
            localsystem = system.filter_by_cell(cell) # Get walkers in the cell
            weights    = localsystem.weights
            walkers    = localsystem.walkers

            print time.asctime(), 'Resampling cell', cell, len(walkers) ,'walkers'

            # Only resample if there exist walkers in the cell
            if not len(walkers) > 0: continue

            ### sort the walkers in descending order based on their weights,
            ##+ this ensures only walkers whose weight > targetWeight are split.
            # This method returns a list of indices in the supplied list/array
            # that correspond to the sorted order of the original.
            mywalkers = list(argsort(-weights))
            print '\tmywalkers:', mywalkers

            ### sanity check
            # Ensure that the list of weights was sorted correctly (i.e., the
            # largest weight comes first)
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
            # The above "continue" check ensures that there is at least one
            # available to the algorithm. Remember that pop takes from the
            # back of a list, so we start with the smallest weight.
            x = mywalkers.pop()

            ### keep track of active walkers for splitting
            activewalkers = 0

            ### The algorithm terminates since the last walker is removed
            ##+ from 'mywalkers' when W == tw. The max number of
            ##+ iterations is bounded by
            ##+ 'len(where(system.cell == cell) + targetwalkers'
            while True: # exit using break

                Wx = weights[x]
                currentWalker = walkers[x]
                print '\tweight of', x, 'is', Wx

                # Split walkers with weight geq the target weight
                # and always split the last walker.
                if Wx >= tw or len(mywalkers) == 0:

                    ### choose number of times to split
                    ##+ r = floor( Wx / tw )
                    ### min, max: work around round-off errors
                    # Per the algorithm, a walker is split r times
                    # corresponding to the nearest integer multiple to the
                    # target weight (e.g., if Wx = 3.5 * tw, split 3 times).
                    r = max(1, int(floor( Wx/tw )) )
                    r = min(r, self.targetwalkers - activewalkers)
                    activewalkers += r
                    print '\tactive walkers', activewalkers

                    ### split the current walker
                    # Note: if r <= 0, repeat will occur 0 times.
                    # This takes advantage of a strange bug in the python
                    # source that sets negative times values equal to 0 in
                    # itertools.repeat when times is not set by keyword.
                    # CHANGE CHANGE CHANGE CHANGE CHANGE
                    print '\tsplitting', x, r, 'times'
                    for _ in itertools.repeat(x, r):
                        # Add a new walker with the target weight to the system
                        # for each time that the walker must be split. Update
                        # the output file as well.
                        w = currentWalker.restart(weight=tw)
                        newsystem.add_walker(w)
                        histfile_fd.write(str(w.initid)+','+str(currentWalker.id)+','+str(w.id)+'\n')


                    ### update the weights for the current walker and mark
                    ##+ for reconsideration
                    # Decrease the current weight to be less than or equal to
                    # the target weight if there are not enough walkers for the
                    # new system and put the walker back on the evaluation list
                    if activewalkers < self.targetwalkers and Wx - r * tw > 0.0:
                        mywalkers.append(x)
                        weights[x] = Wx - r * tw
                        print '\tupdated weights of', x

                    ### continue the loop?
                    # Keep looping until all walkers have been processed.
                    if len(mywalkers) > 0:
                        x = mywalkers.pop()
                    else: break

                ### merge
                # Merge walkers if the weight is less than the target.
                else:
                    # Get another walker to merge with
                    y = mywalkers.pop()
                    print '\tmerging', x, y

                    # Set the merged weight to be the sum of the two weights.
                    Wy = weights[y]
                    Wxy = Wx + Wy

                    # Randomly choose a walker to continue evaluating. This
                    # effectively removes one of them from the list (i.e. one
                    # was absorbed).
                    p = np.random.random()
                    if p < Wy / Wxy:
                        x = y
                    
                    # Reset the weight of the walker under evaluation to the
                    # combined weight.
                    weights[x] = Wxy

                # NOTE: if the walker was merged, a new one is not popped for
                # evaluation in the next loop iteration. Instead, keep working
                # with whichever was kept at the end of the merge section.
        # Once all walkers have been processed, return the new system, which
        # should have all walkers initialized with the target weight.
        histfilefd.close()
        return newsystem

class MultiColor(OneColor):

    """
    awe.resample.MultiColor
    """

    def __init__(self, nwalkers, partition):
        OneColor.__init__(self, nwalkers)
        self.partition   = partition
        ncolors          = partition.ncolors
        self.transitions = np.zeros((ncolors, ncolors))
	self.iteration = 1

        self.cellweights_path = os.path.join(OUTPUT_DIR, 'cell-weights.csv')
        makedirs_parent(self.cellweights_path)
	of = open(self.cellweights_path,'a')
	of.write('%iteration,cellid,color,total_weight \n')
        of.close()

        self.tmat_path = os.path.join(OUTPUT_DIR, 'color-transition-matrix.csv')
        self.tmat_header = textwrap.dedent(
            """\
            # An ((N+1)*2, 2, 2) matrix, where N is the number of iterations
            # Can be loaded like:
            #   m = np.loadtxt("color-transitions-matrix.csv",delimiter=",")
            #   itrs = m.shape[0] / 2
            #   T = m.reshape((itrs, 2, 2))
            """
            )


    def resample(self, system):

        ### update colors
        ncolors          = self.partition.ncolors
        trans = np.zeros((ncolors,ncolors))

        for w in system.walkers:
            cell     = system.cell(w.assignment)
            oldcolor = w.color
            newcolor = cell.core

            # sanity check: all walkers must have a color, but not all cells have a core.
            assert w.color is not None
            assert w.color >= 0

            if not cell.core == aweclasses.DEFAULT_CORE and not w.color == cell.core:
                oldcolor = w.color
                newcolor = cell.core
                print 'Updating color:', w, oldcolor, '->', newcolor
                w.color = newcolor
            else:
                oldcolor = newcolor = w.color

            trans[oldcolor, newcolor] += w.weight
        self.transitions = np.append(self.transitions,trans,axis=0)

        ### resample individual colors using OneColor algorithm
        newsystem = aweclasses.System(topology=system.topology)
        for color in system.colors:
            thiscolor  = system.filter_by_color(color)
            print time.asctime(), 'Resampling color', color, len(thiscolor.walkers), 'walkers'
            resampled  = OneColor.resample(self, thiscolor)
            newsystem += resampled

        of = open(self.cellweights_path,'a')
        for cell in newsystem.cells:
	    thiscell = system.filter_by_cell(cell)
	    for color in thiscell.colors:
	        thiscolor = thiscell.filter_by_color(color)
		of.write(str(self.iteration)+','+str(cell.id)+','+str(color)+','+str(sum(thiscolor.weights))+'\n')
    	of.close()
	    self.iteration += 1

        self.save_transitions(self.tmat_path)

        return newsystem

    def save_transitions(self, path):
        print time.asctime(), 'Saving transition matrix to', repr(path)
        fd = open(path, 'w')
        try:
            fd.write(self.tmat_header)
            np.savetxt(fd, self.transitions, delimiter=',')
        finally:
            fd.close()

class SuperCell(MultiColor):
    def __init__(self,nwalkers,partition,cellmapf):
        MultiColor.__init__(self,nwalkers,partition)
        self.cellmap = []
        for line in open(cellmapf):
            self.cellmap.append(int(line))

    def resample(self,system):
        cellmap = self.cellmap
        for w in system.walkers:
            w.assignment = cellmap[w.assignment]
        newsystem = MultiColor.resample(self,system)
        tmat = os.path.join(OUTPUT_DIR, 'transition-matrix.csv')
        makedirs_parent(tmat)
        MultiColor.save_transitions(tmat)
        return newsystem

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
    def __init__(self, resampler, datfile='isaver.dat'):
        self.resampler = resampler
        self.datfile   = datfile
        self.iteration = 0

    @typecheck(aweclasses.System, mode=str)
    def save(self, system, mode='a'):
        if self.iteration == 0:
            with open(self.datfile, 'a') as fd:
                fd.write(self.heading())
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

        newsystem = self.resampler.resample(system)
        self.iteration += 1
        self.save(newsystem, mode='a')

        return newsystem

class SaveWeights(ISaver):

    def __init__(self, resampler, datfile=None):
        datfile = datfile or os.path.join(OUTPUT_DIR, 'walker-weights.csv')
        ISaver.__init__(self, resampler, datfile=datfile)

    def _heading(self):
        return \
            '# Each line represents a walker at:\n' + \
            '# walkerid,iteration,cell,weight,color\n'

    def _save(self, system, mode='a'):
        print time.asctime(), 'Saving weights to', self.datfile

        ### all the walkers in a cell have the same weight, so we only
        ### need to save the walkerid, iteration, cell, weight, and color for each walker

        with open(self.datfile, mode) as fd:
            for i, w in enumerate(system.walkers):
                s = '%(wid)d,%(iteration)d,%(cell)d,%(weight)f,%(color)s\n' % {
                    'wid'       : w.id           ,
                    'iteration' : self.iteration ,
                    'cell'      : w.assignment   ,
                    'weight'    : w.weight       ,
                    'color'     : w.color } # if w.color is not None else 'nan'       }
                fd.write(s)
