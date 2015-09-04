#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Voronoi diagram from a list of points
# Copyright (C) 2011  Nicolas P. Rougier
#
# Distributed under the terms of the BSD License.
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def circumcircle(P1, P2, P3):
    '''
    Return center of the circle containing P1, P2 and P3

    If P1, P2 and P3 are colinear, return None

    Adapted from:
    http://local.wasp.uwa.edu.au/~pbourke/geometry/circlefrom3/Circle.cpp
    '''
    delta_a = P2 - P1
    delta_b = P3 - P2
    if np.abs(delta_a[0]) <= 0.000000001 and np.abs(delta_b[1]) <= 0.000000001:
        center_x = 0.5*(P2[0] + P3[0])
        center_y = 0.5*(P1[1] + P2[1])
    else:
        aSlope = delta_a[1]/delta_a[0]
        bSlope = delta_b[1]/delta_b[0]
        if np.abs(aSlope-bSlope) <= 0.000000001:
            return None
        center_x= (aSlope*bSlope*(P1[1] - P3[1]) + bSlope*(P1[0] + P2 [0]) \
                        - aSlope*(P2[0]+P3[0]))/(2.*(bSlope-aSlope))
        center_y = -(center_x - (P1[0]+P2[0])/2.)/aSlope + (P1[1]+P2[1])/2.
    return center_x, center_y

def voronoi(X,Y):
    ''' Return line segments describing the voronoi diagram of X and Y '''
    P = np.zeros((X.size+4,2))
    P[:X.size,0], P[:Y.size,1] = X, Y
    # We add four points at (pseudo) "infinity"
    m = max(np.abs(X).max(), np.abs(Y).max())*1e5
    P[X.size:,0] = -m, -m, +m, +m
    P[Y.size:,1] = -m, +m, -m, +m
    D = matplotlib.tri.Triangulation(P[:,0],P[:,1])
    T = D.triangles

    #axes = plt.subplot(1,1,1)
    #plt.scatter(X,Y, s=5)
    #patches = []
    #for i,triang in enumerate(T):
    #    polygon = Polygon(np.array([P[triang[0]],P[triang[1]],P[triang[2]]]))
    #	patches.append(polygon)
    #colors = 100*np.random.rand(len(patches))
    #p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    #p.set_array(np.array(colors))
    #axes.add_collection(p)
    #plt.colorbar(p)
    ##lines = matplotlib.collections.LineCollection(segments, color='0.75')
    ##axes.add_collection(lines)
    #plt.axis([0,1,0,1])
    #plt.show()
    #plt.savefig('test0.png')

    n = T.shape[0]
    C = np.zeros((n,2))
    for i in range(n):
        C[i] = circumcircle(P[T[i,0]],P[T[i,1]],P[T[i,2]])
    X,Y = C[:,0], C[:,1]
    #segments = []
    #for i in range(n):
    #    for k in D.neighbors[i]:
    #        if k != -1:
    #            segments.append([(X[i],Y[i]), (X[k],Y[k])])
    cells = [[] for i in range(np.max(T)+1)]
    for i,triang in enumerate(T):
        cells[triang[0]].append([X[i],Y[i]])
        cells[triang[1]].append([X[i],Y[i]])
	cells[triang[2]].append([X[i],Y[i]])
    for i,cell in enumerate(cells):
        angle = []
        for coord in cell:
            angle.append(np.arctan2((coord[1]-P[i,1]),(coord[0]-P[i,0])))
        id = np.argsort(-np.array(angle))
        cells[i] = np.array([cell[j] for j in id])
    return cells

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    P = np.random.random((2,256))
    X,Y = P[0],P[1]
    fig = plt.figure(figsize=(10,10))
    axes = plt.subplot(1,1,1)
    plt.scatter(X,Y, s=5)
    #segments = voronoi(X,Y)
    cells = voronoi(X,Y)
    patches = []
    for cell in cells:
        polygon = Polygon(cell,True)
        patches.append(polygon)
        #plt.scatter(cell[:,0],cell[:,1])
    colors = 100*np.random.rand(len(patches))
    print(colors)
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    p.set_array(np.array(colors))
    axes.add_collection(p)
    plt.colorbar(p)
    #lines = matplotlib.collections.LineCollection(segments, color='0.75')
    #axes.add_collection(lines)
    plt.axis([0,1,0,1])
    plt.show()
    plt.savefig('test.png') 
