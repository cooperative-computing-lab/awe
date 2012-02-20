#!/usr/bin/env python

from sys import argv
import sys
import os
import os.path
import shutil

OUTPUT_DIR = "output"
SCRATCH_DIR = "scratch"


if __name__ == "__main__":
    if len(argv) != 4:
        print "Usage: %s input.pdb topol.top sim.mdp" % argv[0]
        sys.exit(1)
    
    inputStruct = argv[1]
    topol = argv[2]
    simMdp = argv[3]
    
    inBaseName = os.path.basename(inputStruct)
    inBasePrefix = inBaseName[:inBaseName.rfind(".")]
    scratchBase = os.path.join(SCRATCH_DIR, inBasePrefix)

    # change NLYP to NLYS
    # NLYP was used in the old FFAmber -- NLYS is used now
    updStruct = scratchBase + ".fixed.pdb"
    if os.path.exists(updStruct):
        os.remove(updStruct)
    #cmd = "cat %(inputStruct)s | sed -e 's/LYS /NLYS/g' | sed -e 's/NLYP/NLYP/g' > %(updStruct)s" % locals()
    shutil.copy(inputStruct, updStruct)
    #print cmd
    #os.system(cmd)
    
    # grompp for sim
    simOutMdp = scratchBase + ".sim.out.mdp"
    simTpr = scratchBase + ".tpr"
    cmd = "grompp -f %(simMdp)s -po %(simOutMdp)s -c %(updStruct)s -p %(topol)s -o %(simTpr)s" % locals()
    if os.system(cmd) != 0:
        sys.exit(1)

    shutil.copy(simTpr, OUTPUT_DIR)
    
