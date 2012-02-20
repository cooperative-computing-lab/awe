
from awesome   import AWE, Cells, Topology
from somewhere import loadcoords

import glob


### coordinates need a topology
topol      = Topology.from_pdb('atoms.pdb')

### load cells, in this case from an msm
cellcoords = load_msm_states('states.hd5')
cells      = Cells(cellcoords)


### associate cells with colors
##+ this example is simplistic, but could conceivably
##+ be done via a more complex procedure
##+ (ie RMSD of cell center to some reference)
cells.color('unfolded', [0,1,2])
cells.color('folded'  , [3,4,5])
cells.color('unfolded', [6,7])
cells.color('folded',   [8,9])


### associate initial walkers to the cells
##+ one walker per cell, but this allows arbitrary association
for cellid, pdb in enumerate(glob.glob('walkers/*.pdb')):
    coords = loadcoords(pdb)
    cells.walker(cellid, coords)


### setup and run the AWE method
awe = AWE(cells, topol=topol, wq=wq, iterations=42, name='awesome!')

# worker script
awe.execute('md_and_assign.sh')

# cache these files on the worker
awe.cache('states.hd5')
awe.cache('Assign.py')

# save the state after each call to 'resample'
awe.logto('log.dat')

# run several iterations of the MD+Assign, Resample procedure
awe.run()
