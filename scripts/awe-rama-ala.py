#!/usr/bin/env python

from prody import *
import os
import sys
import numpy as np
import awe.voronoi as vn
from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

setVerbosity('none')
#input: pdb file name, cell.dat file name
#output: phi and psi angels in a two dimensional list
def parseCenter(pdbf,cellf):
    ala = parsePDB(pdbf)
    phi = []
    psi = []
    coords = []
    res = []
    for line in open(cellf,'r'):
        try:
	    coord = float(line)
	    res.append(coord)
	except ValueError:
	    continue
	if len(res) == 3:
	    coords.append(res)
	    res = []
	if len(coords) == 22:
            ala.setCoords(np.array(coords))
	    writePDB('tmp.pdb',ala)
	    os.system('g_chi -s tmp.pdb -f tmp.pdb -rama >& tmp.log')
            phipsif = open('ramaPhiPsiALA2.xvg')
	    fcontent = phipsif.readlines()
	    phipsil = fcontent[len(fcontent)-1]
	    phipsitext = phipsil.split()
	    phi.append(float(phipsitext[0]))
            psi.append(float(phipsitext[1]))
            phipsif.close()
            os.system('rm tmp.pdb chi.log ramaPhiPsiALA2.xvg tmp.log order.xvg')
	    #phipsi.append([calcPhi(alacenter),calcPsi(alacenter)])
	    coords = []
    return np.array([phi,psi])

#phipsi = parseCenter('testinput/state0.pdb','testinput/cells.dat')

def extractWeight(weightf,ncell):
    weights = np.zeros([ncell])
    l = 1;
    for line in open(weightf,'r'):
        if l > 1:
	    content = line.split(',')
            cellid = int(content[1])
	    w = float(content[3])
	    weights[cellid] += w
        l += 1
    return weights
    
def ramaColor(weights,maxfe):
    weights = weights/sum(weights)
    for i in range(len(weights)):
        if weights[i] != 0:
            weights[i] = (-1)*np.log(weights[i])
    maxb = np.max(weights)
    if maxb > maxfe:
        maxb = maxfe
    for i in range(len(weights)):
        if weights[i] == 0:
	    weights[i] = maxb
    return weights

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-w','--weigthFile',dest='wf',help='file with weight of cells',default='output.dat')
    #parser.add_option('-m','--maxfe',dest='maxfe',help='maxfe for Ramachandran plot',type='float',default=1000)
    parser.add_option('-p','--pdbfile',dest='pdbf',help='pdb file of protein structure')
    parser.add_option('-c','--cellsCenterFile',dest='cellf',help='file of coordinate of cells center',default='testinput/cells.dat')
    parser.add_option('-n','--numberCell',dest='ncell',help='number of cells in AWE',type='int')

    args = parser.parse_args()
    weightf = args[0].wf
    #maxfe = args[0].maxfe
    pdbf = args[0].pdbf
    cellf = args[0].cellf
    ncell = args[0].ncell

    maxfe = 10
    print 'Converting coordinates to phi-psi angles'
    cellcenter = parseCenter(pdbf,cellf)
    print 'Reading cell weights'
    cellw = extractWeight(weightf,ncell)
    colorw = ramaColor(cellw,maxfe)

    print 'Plotting Ramachandran'
    X = cellcenter[0]
    Y = cellcenter[1]
    cells = vn.voronoi(X,Y)
    patches = []
    for cell in cells:
        polygon = Polygon(cell,True)
        patches.append(polygon)

    p = PatchCollection(patches,cmap=matplotlib.cm.jet,alpha=0.4)
    p.set_array(colorw)

    fig = plt.figure(figsize=(10,10))
    axes = plt.subplot(1,1,1)
    plt.scatter(X,Y)
    axes.add_collection(p)
    plt.colorbar(p)
    plt.savefig('awe-rama-ala.png')


