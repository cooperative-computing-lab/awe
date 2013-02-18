# -*- mode: Python; indent-tabs-mode: nil -*-  #

import awe

import numpy as np
import sys
import getopt
import os
import getpass

#-----Simulation Default Values-----
iterations = 5
nwalkers   = 4
nstates    = 100
restarts   = float('inf')
maxreps    = 50

#-----WQ Default Values-----
wq_port = 0
wq_fast_abort_multiplier = -1.0
wq_proj_name = None 
wq_debug_flag = None

#-----Main Program------
if __name__ == "__main__":

	help_str = "Usage: %s\n" % sys.argv[0]
	help_str += "-i		-	specify the number of iterations (default: %d)\n" % iterations
	help_str += "-w		-	specify the number of walkers (default: %d)\n" % nwalkers
	help_str += "-s		-	specify the number of states (default: %d)\n" % nstates
	help_str += "-r		-	specify the number of restarts allowed (default: %f)\n" % restarts
	help_str += "-p		-	specify the port to use for Work Queue (default: arbitrary)\n"
	help_str += "-n		-	specify a project name for Work Queue to use\n"
	help_str += "-f		-	specify the Work Queue fast abort multipler\n"
	help_str += "-d		-	print Work Queue debug messages \n"
	help_str += "-h		-	help"

	try:
		opts, args = getopt.getopt(sys.argv[1:], "d:f:hi:n:p:r:s:w", ["help"])
	except getopt.GetoptError, err:
		print str(err) 
		print help_str
		sys.exit(1)

	#Parse command line arguments.
	for o, a in opts:
		if o in ("-d"):
			wq_debug_flag = a 
		elif o in ("-f"):
			wq_fast_abort_multiplier = float(a) 
		elif o in ("-h", "--help"):
			print help_str
			sys.exit(0)
		elif o == "-i":
			iterations = int(a)
		elif o in ("-n"):
			wq_proj_name = a
		elif o in ("-p"):
			wq_port = int(a) 
		elif o in ("-r"):
			restarts = int(a) 
		elif o in ("-s"):
			nstates = int(a) 
		elif o in ("-w"):
			nwalkers = int(a)

	cfg           = awe.workqueue.Config()
	cfg.fastabort = wq_fast_abort_multiplier 
	cfg.restarts  = restarts 
	cfg.maxreps   = maxreps
	cfg.name      = wq_proj_name 
	cfg.port      = wq_port
	
	if wq_debug_flag:
		cfg.debug     =	wq_debug_flag 

        # The "main" function of the worker
	cfg.execute('testinput/execute-task.sh')

        # Binaries to run MD and assignment steps
	cfg.cache('awesetup/binaries/$OS-$ARCH/pdb2gmx')
	cfg.cache('awesetup/binaries/$OS-$ARCH/grompp')
	cfg.cache('awesetup/binaries/$OS-$ARCH/mdrun')
	cfg.cache('awesetup/binaries/$OS-$ARCH/assign')

	cfg.cache('awesetup/gmxtopologies')  # required for running gromacs for MD
	cfg.cache('testinput/sim.mdp') # Gromacs simulation parameters
	cfg.cache('testinput/env.sh') # setting up the worker execution environment
	cfg.cache('testinput/cells.dat') # cell definitions
	cfg.cache('testinput/CellIndices.dat') # cell atoms to use when assigning
	cfg.cache('testinput/StructureIndices.dat') # walker atoms to use when assigning

        # initialize the weights randomly
	weights   = np.random.random((nstates,nwalkers))
	weights  /= np.sum(weights.flatten())

        # load a topology file
	system    = awe.System(topology = awe.PDB('testinput/state0.pdb'))

        # 2-color awe needs states assigned to a region
	partition = awe.SinkStates()
	partition.add(0, *range(0,nstates/2))
	partition.add(1, *range(nstates/2,nstates))

        # load the initial cells and walkers
	srcdir = 'awesetup/pdbs/ala'
	for i in xrange(nstates):

	    if i < nstates / 3:
		cell = awe.Cell(i, core=0)
	    elif i > 2 * nstates / 3:
		cell = awe.Cell(i, core=1)
	    else:
		cell = awe.Cell(i)

	    color = partition.color(cell)
	    system.add_cell(cell)


	    for j in xrange(nwalkers):

		pdbpath = os.path.join(srcdir, 'State%d-%d.pdb' % (i, j))
		pdb     = awe.PDB(pdbpath)
		w       = awe.Walker(start=pdb.coords, assignment=i, color=color, weight=weights[i,j], cellid=cell.id)
		system.add_walker(w)

        # define the AWE resampling algorithm to use
	multicolor = awe.resample.MultiColor(nwalkers, partition)
	resample   = awe.resample.SaveWeights(multicolor)
	adaptive   = awe.AWE( wqconfig   = cfg,
			      system     = system,
			      iterations = iterations,
			      resample   = resample,
			      checkpointfreq = 1)

	adaptive.run()

	print 'Run time:', awe.time.time(), 's'
	sys.exit(0)
