# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""

###############################################################################
# NOTE - Jeff Kinnison, July 2015
#
# This file was apparently written using text editors with different tab
# lengths and space/tab schemes. I have tried to fix this, but it may still
# require some work before it will run. There isn't a great way to do this
# other than manually thanks to Python's scoping by indentation.
#
# Run the file, and it will tell you about any overt errors. Bugs may still
# propagate from indentation where there should be none, however. I used
# context to determine how tabs should be fixed but may have gotten it wrong
# in places.
# 
# YOU HAVE BEEN WARNED
#
# 5102 yluJ ,nosinniK ffeJ - ETON
###############################################################################

from .util import typecheck, returns, makedirs_parent
from . import aweclasses

import numpy as np
import itertools
import time, os
import textwrap

OUTPUT_DIR = os.getcwd()


class IResampler(object):
    """
    Interface that all resampling methods should implement.

    Fields:
        None

    Methods:
        resample - divide and kill off walkers until weights are equal to
                   the group's mean weight
    """

    def resample(self, walkers):
        """
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
        print(time.asctime(), 'Resampling')
        return self.resample(s1)


class Identity(IResampler):
    """
    A resampler that does not alter the System it receives.

    Fields:
        None

    Methods:
        resample - returns the identity of the System
    """

    def resample(self, walkers):
        """
        Perform the identity operation on the System.

        Parameters:
            walkers - an awe.aweclasses.System instance

        Returns:
            The supplied awe.aweclasses.System instance
        """

        return walkers


class OneColor(IResampler):
    """
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

            print(time.asctime(), 'Resampling cell', cell, len(walkers) ,'walkers')

            # Only resample if there exist walkers in the cell
            if not len(walkers) > 0: continue

            ### sort the walkers in descending order based on their weights,
            ##+ this ensures only walkers whose weight > targetWeight are split.
            # This method returns a list of indices in the supplied list/array
            # that correspond to the sorted order of the original.
            mywalkers = list(argsort(-weights))
            print('\tmywalkers:', mywalkers)

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
            print('\tW', W, 'tw', tw)

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
                print('\tweight of', x, 'is', Wx)

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
                    print('\tactive walkers', activewalkers)

                    ### split the current walker
                    # Note: if r <= 0, repeat will occur 0 times.
                    # This takes advantage of a strange bug in the python
                    # source that sets negative times values equal to 0 in
                    # itertools.repeat when times is not set by keyword.
                    # CHANGE CHANGE CHANGE CHANGE CHANGE
                    print('\tsplitting', x, r, 'times')
                    for _ in itertools.repeat(x, r):
                        # Add a new walker with the target weight to the system
                        # for each time that the walker must be split. Update
                        # the output file as well.
                        w = currentWalker.restart(weight=tw)
                        newsystem.add_walker(w)
                        histfile_fd.write(str(w.initid)+','+ \
                            str(currentWalker.id)+','+str(w.id)+'\n')


                    ### update the weights for the current walker and mark
                    ##+ for reconsideration
                    # Decrease the current weight to be less than or equal to
                    # the target weight if there are not enough walkers for the
                    # new system and put the walker back on the evaluation list
                    if activewalkers < self.targetwalkers and Wx - r * tw > 0.0:
                        mywalkers.append(x)
                        weights[x] = Wx - r * tw
                        print('\tupdated weights of', x)

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
                    print('\tmerging', x, y)

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
    A resampler for handling multiple macro-states (metastable states)
    by following transitions between states and resampling each
    individually.

    Fields:
        targetwalkers    - the number of walkers, creates a target weight
        histfile         - the path to walker-history.csv
        partition        - awe.aweclasses.SinkStates instance
        transitions      - a matrix containing transition rates between
                           macro-states
        iteration        - the current iteration of the 
        cellweights_path - the path to cell-weights.csv
        tmat_path        - the path to color-transition-matrix.csv
        tmat_header      - the header for the color-transition-matrix.csv file


    Methods:
        resample         - resample across multiple states
        save_transitions - output the state transition matrix
    """

    def __init__(self, nwalkers, partition):
        """
        Initialize a new instance of MultiColor.

        Parameters:
            nwalkers  - the target number of walkers
            partition - an aweclasses.SinkStates instance

        Returns: 
            None
        """

        # Set up the output filepath and target weights
        OneColor.__init__(self, nwalkers)

        # Partitioning of the states (e.g. conformation space cells)
        self.partition   = partition
        ncolors          = partition.ncolors
        
        # Matrix to store transitions between states
        self.transitions = np.zeros((ncolors, ncolors))
	    self.iteration = 1

        # PAths to output files and necessary header information
        self.cellweights_path = os.path.join(OUTPUT_DIR, 'cell-weights.csv')
        makedirs_parent(self.cellweights_path)
    	of = open(self.cellweights_path,'a')
    	of.write('%iteration,cellid,color,total_weight \n')
        of.close()

        self.tmat_path = os.path.join(OUTPUT_DIR,
                'color-transition-matrix.csv')
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
        """
        Resample over multiple cells by calling the OneColor algorithm on each
        in succession. Note: look into parallelizing running the OneColor algo
        on each cell.

        Parameters:
            system - a processed instance of aweclasses.System to be resampled

        Returns:
            A new aweclasses.System instance representing the resampled input
            system
        """

        ### update colors
        # Create a new transition matrix representing the current set of
        # macro-states.
        ncolors = self.partition.ncolors
        trans = np.zeros((ncolors,ncolors))

        # Process each walker in the system
        for w in system.walkers:
            # Determine the previous and current states of the walker
            cell     = system.cell(w.assignment)
            oldcolor = w.color
            newcolor = cell.core

            # sanity check: all walkers must have a color, but not all cells have a core.
            # Make sure that all walkers are in a state
            assert w.color is not None
            assert w.color >= 0

            # Update the walker to reflect its current state
            # Most of this seems redundant from above.
            if not cell.core == aweclasses.DEFAULT_CORE and not w.color == cell.core:
                oldcolor = w.color   # This line seems redundant
                newcolor = cell.core # This line seems redundant
                print('Updating color:', w, oldcolor, '->', newcolor)
                w.color = newcolor # Here we actually update the walker
            else:
                # This is redundant. If the two are already equal the
                # assignment changes nothing, and this branch is only accessed
                # when the two are equal.
                oldcolor = newcolor = w.color

            # Add the transition to the transition matrix.
            trans[oldcolor, newcolor] += w.weight
        
        # Add all transitions to the instance transition matrix
        self.transitions = np.append(self.transitions,trans,axis=0)

        ### resample individual colors using OneColor algorithm
        newsystem = aweclasses.System(topology=system.topology)
        
        # Do resampling per state using the above OneColor algorithm
        for color in system.colors:
            # Get all walkers of the color
            thiscolor  = system.filter_by_color(color)
            print(time.asctime(), 'Resampling color', color, len(thiscolor.walkers), 'walkers')
            
            # Perform resampling and add to the new system
            resampled  = OneColor.resample(self, thiscolor)
            newsystem += resampled

        # Output information about the new system for restart and analysis use
        of = open(self.cellweights_path,'a')
        for cell in newsystem.cells:
            # Get the corresponding cell from the old system
    	    thiscell = system.filter_by_cell(cell)
    	    
            for color in thiscell.colors:
                # Get the corresponding color from the cell in the old system
    	        thiscolor = thiscell.filter_by_color(color)

                # Output information about the cell and color
    	        of.write(str(self.iteration)+','+str(cell.id)+','+str(color)+ \
                    ','+str(sum(thiscolor.weights))+'\n')
        of.close()

        # This iteration is finished, so increment
        self.iteration += 1

        # Output the transition matrix
        self.save_transitions(self.tmat_path)

        return newsystem

    def save_transitions(self, path):
        """
        Output the state transition matrix to a file.

        Parameters:
            None

        Returns:
            None
        """

        print(time.asctime(), 'Saving transition matrix to', repr(path))
        fd = open(path, 'w')
        
        try:
            # Output the transitions from this iteration only
            fd.write(self.tmat_header)

            # This saves the current transitions in case of failure
            np.savetxt(fd, self.transitions, delimiter=',')
        
        finally:
            fd.close()




############# END OF CLASSES USED BY BADI's AWE-WQ IMPLEMENTATION #############

# Classes beyond this point were used in the release version od AWE and are
# either legacy classes or those to be implemented in the future.



class SuperCell(MultiColor):
    """
    THIS DOES NOT SEEM TO BE USED. There is not much context for determining
    the purpose of this class, so any explanations for additions to the
    MultiColor class are interpretations of the code.

    This class seems to allow saving the state partitioning scheme into an
    arbitrary list, and then walkers are assigned cells from that list. It acts
    as an abstraction for determining which cells walkers are in, however that
    abstraction seems to break down after one iteration. Possibly intended to
    restart AWE from the last iteration it completed by re-performing
    resampling based on a saved cell state.

    Fields: 
        targetwalkers    - the number of walkers, creates a target weight
        histfile         - the path to walker-history.csv
        partition        - awe.aweclasses.SinkStates instance
        transitions      - a matrix containing transition rates between
                           macro-states
        iteration        - the current iteration of the 
        cellweights_path - the path to cell-weights.csv
        tmat_path        - the path to color-transition-matrix.csv
        tmat_header      - the header for the color-transition-matrix.csv file
        cellmap          - 

    Methods:
        resample         - 
        save_transitions -
    """

    def __init__(self,nwalkers,partition,cellmapf):
        """
        Initialize a new instance of SuperCell.

        Parameters:
            nwalkers  - the target number of walkers for the system
            partition - a partitioning of the system states
            cellmapf  - a file containing the mapping of walkers into cells?
        """

        MultiColor.__init__(self,nwalkers,partition)
        self.cellmap = []
        for line in open(cellmapf):
            self.cellmap.append(int(line))

    def resample(self,system):
        """
        A resampling that pulls walker cell assignments from an external file
        and then resamples the whole system under MultiColor. Given that the
        assignments are left in post-mapping state at the end of the method,
        it cannot be used multiple times.

        Parameters:
            system - a processed aweclasses.System object

        Returns:
            A new aweclasses.System instance representing the resampled input 
            system
        """
        
        cellmap = self.cellmap
        
        # Determine the cell that each walker is assigned to so that
        # cells can be stored in some arbitrary order (?)
        for w in system.walkers:
            w.assignment = cellmap[w.assignment]

        # Perform MultiColor resampling
        newsystem = MultiColor.resample(self,system)
        
        # Save the transition matrix
        tmat = os.path.join(OUTPUT_DIR, 'transition-matrix.csv')
        makedirs_parent(tmat)
        MultiColor.save_transitions(tmat)

        return newsystem

class IPlotter(IResampler):
    """
    NOTE: This class uses legacy nomenclature "walkergroup", so it was likely
    never actually used. However, a utility for plotting attributes of a system
    is definitely interesting.

    The interface for a plotting utility with resampling capability.

    Fields:
        plotfile - the file to which the plot should be saved

    Methods:
        resample - potential method for resampling
        compute  - potential function for computing plot variables
        plot     - potential method for plotting computed variables
    """

    def __init__(self, **kws):
        """
        Compute some attribute of the system passed in.

        Parameters:
            kws - keyword arguments
                  "plotfile": a file to save the plot to

        Returns:
            None
        """

        self.plotfile = kws.pop('plotfile', 'plot.png')

    def compute(self, walkergroup):
        """
        Compute some attribute of the system passed in.

        Parameters:
            walkergroup - a system for which attributes should be computed

        Raises:
            NotImplementedError

        Returns:
            None
        """

        raise NotImplementedError

    def plot(self):
        """
        Plot some computed attribute of the system passed in.

        Parameters:
            None

        Raises:
            NotImplementedError

        Returns:
            None
        """

        raise NotImplementedError

    def __call__(self, walkers):
        ws = IResampler.__call__(self, walkers)
        self.compute(ws)
        self.plot()
        return ws

class ISaver(IResampler):
    """
    Utility for saving information after each resampling.

    Fields:
        resampler - a subclass of IResampler
        datfile   - the filepath to which data should be saved
        iteration - the iteration of the resampling algorithm

    Methods:
        resample - resample a system by the resampler
        save     - save data to datfile
        heading  - the heading for the datfile
    """

    @typecheck(IResampler, datfile=str)
    def __init__(self, resampler, datfile='isaver.dat'):
        """
        Initialize a new instance of ISaver.

        Parameters:
            resampler - a resampler class to use for the actual resampling
            datfile   - the file to which data should be saved

        Returns:
            None
        """

        self.resampler = resampler
        self.datfile   = datfile
        self.iteration = 0

    @typecheck(aweclasses.System, mode=str)
    def save(self, system, mode='a'):
        """
        Save data to datfile.

        Parameters:
            system - an instance of aweclasses.System to be saved
            mode   - the mode in which to open datfile
        """

        if self.iteration == 0:
            with open(self.datfile, 'a') as fd:
                fd.write(self.heading())
        self._save(system, mode=mode)

    @returns(str)
    def heading(self):
        """
        Get the heading for the output data.

        Parameters:
            None

        Returns:
            The heading for the output file
        """

        return self._heading()

    def _save(self, system, mode='a'):
        raise NotImplementedError

    def _heading(self):
        return ''

    @typecheck(aweclasses.System)
    @returns(aweclasses.System)
    def resample(self, system):
        """
        Call the resampler and save the new system to file.

        Parameters:
            system - an instance of aweclasses.System to be resampled and saved

        Returns:
            A new aweclasses.System instance representing the resampled input 
            system
        """

        newsystem = self.resampler.resample(system)
        self.iteration += 1
        self.save(newsystem, mode='a')

        return newsystem

class SaveWeights(ISaver):
    """
    A utility for saving walker weights to a file.

    Fields:
        resampler - a subclass of IResampler
        datfile   - the filepath to which data should be saved
        iteration - the iteration of the resampling algorithm

    Methods:
        resample - resample a system by the resampler
        save     - save data to datfile
        heading  - the heading for the datfile
    """

    def __init__(self, resampler, datfile=None):
        """
        Initialize a new instance of SaveWeights.

        Parameters:
            resampler - a subclass of IResampler
            datfile   - the filepath to which data should be saved
        """

        datfile = datfile or os.path.join(OUTPUT_DIR, 'walker-weights.csv')
        ISaver.__init__(self, resampler, datfile=datfile)

    def _heading(self):
        return \
            '# Each line represents a walker at:\n' + \
            '# walkerid,iteration,cell,weight,color\n'

    def _save(self, system, mode='a'):
        print(time.asctime(), 'Saving weights to', self.datfile)

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
