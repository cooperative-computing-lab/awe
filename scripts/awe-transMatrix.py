#!/usr/bin/env python

from optparse import OptionParser
import numpy as np

#build a dict recording walker history from the walkerhistory.dat file
def buildpath(whistory_f):
    whistory_fd = open(whistory_f,'r')
    path = {}
    for line in whistory_fd:
        if line[0] == '%':
            continue
	linecontent = line.split(',')
        try:
	    curnode = linecontent[1]
	    child = int(linecontent[2])
        except IndexError:
            continue
	# the curnode point to a list of its root and children
	try:
            path[curnode].append(child)
        except KeyError:
            path[curnode] = [child]
        path[str(child)] = []
    whistory_fd.close()
    return path

def loadWeightFile(weight_f):
    weight_fd = open(weight_f,'r')
    wcellweight = {}
    for line in weight_fd:
        if line[0] == '#':
	    continue
	else:
	    linecontent = line.split()
	    try:
                wcellweight[linecontent[0]] = [int(linecontent[1]),int(linecontent[2]),float(linecontent[3]),int(linecontent[4])] 
            except IndexError:
                continue
    return wcellweight

#generate transition matrix with a given time lag
def transCount(path,wcellweight,tlag,numcell):
    #tlag defined as number of iterations in the lag time
    #wcellweight defines the walker id, cell, weight and color
    #returns two transition matrices for color 0 and 1
    transC0 = np.zeros([numcell,numcell])
    transC1 = np.zeros([numcell,numcell])
    transC2 = np.zeros([numcell,numcell])
    for walker in path:
        try:
            origcell = wcellweight[walker][1]
            origcolor = wcellweight[walker][3]
        except KeyError:
            continue
        
        children = findChildren(path,[int(walker)],tlag)
	if len(children) == 0:
            continue
        else:
            for w in children:
                strw = str(w)
	        tarcell = int(wcellweight[strw][1])
	        tweight = wcellweight[strw][2]
	        color = wcellweight[strw][3]
                if origcolor == 0:
	            transC0[origcell][tarcell] += tweight
	        else:
	            transC1[origcell][tarcell] += tweight
                #else:
                #    transC2[origcell][tarcell] += tweight
    return (transC0,transC1)
        
#find children of a node with given time lag (apart with given generations)
def findChildren(path,nodes,tlag):
    if tlag <= 0 or len(nodes)==0:
        return nodes
    else:
        newnodes = []
        for w in nodes:
	    try:
		newnodes.extend(path[str(w)])
	    except KeyError:
	        continue
	children=findChildren(path,newnodes,tlag-1)
        return children
	   
def subpath(miniter,maxiter,wcellweight,path):
    newpath = {}
    for w in wcellweight:
        if wcellweight[w][0] >= miniter and wcellweight[w][0] < maxiter:
	    newpath[w] = path[w] 
    return newpath

def buildTrans(assign,tlag,ncell):
    count = np.zeros([ncell,ncell])
    for trj in assign:
        for i in xrange(len(trj)-tlag):
            count[trj[i]][trj[i+tlag]] += 1
    ts = count + np.transpose(count)
    trans = []
    for t in ts:
        if sum(t) != 0:
            trans.append(t/sum(t))
        else:
            trans.append(t)
    return np.array(trans)

def writeCSV(data,filename):
    output = open(filename,'w')
    rows = len(data)
    try:
        cols = len(data[0])
        for i in xrange(rows):
            for j in xrange(cols):
                output.write(str(data[i][j])+',')
            output.write('\n')
    except TypeError:
        for i in xrange(rows):
            output.write(str(data[i])+',')
    output.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-p','--whistoryfile',dest='whistory_f',default='walker-history.csv',help='file that tracks walker id history,defult=walker-history.csv')
    parser.add_option('-w','--weightsfile',dest='wweight_f',default='walker-weights.csv',help='file that records walkerid,iteration,cell,weight,color, defult=walker-weights.csv')
    parser.add_option('-n','--numcell',dest='numcell',type='int',default=-1,help='number of cells')
    parser.add_option('-t','--timelag',dest='tlag',type='int',default=1,help='iterations apart for calculating transition matrix,default=1')
    parser.add_option('-o','--output',dest='ofile',default='./trans.csv',help='output file name, default=./trans.csv')
    parser.add_option('-s','--subpath',dest='splength',type='int',default=-1,help='iteration start from which transition matrix is built,default=0')

    args = parser.parse_args()
    wh_f = args[0].whistory_f
    ww_f = args[0].wweight_f
    ncell = args[0].numcell
    tlag = args[0].tlag
    of = args[0].ofile
    spl = args[0].splength

    print 'Loading walker history file'
    path = buildpath(wh_f)
    print 'Loading walker weights file'
    wcellweight = loadWeightFile(ww_f)

    if spl <= 0:
        print 'Calculate transition matrix'
        trans0,trans1 = transCount(path,wcellweight,tlag,ncell)
	trans = trans0 + trans1

        trans2 = trans+np.transpose(trans)
        transP = []
        for t in trans2:
            if sum(t) != 0:
                transP.append(t/sum(t))
            else:
                transP.append(t)
        writeCSV(transP,of)
        eval,evec = np.linalg.eig(np.array(transP))
        print 'Implied Time Scale (unit iteration length): '+str(-tlag/np.log(eval[1]))
    else:
        print 'Calculate transition matrix'
        newpath = subpath(spl,1000000,wcellweight,path)
        trans0,trans1=transCount(newpath,wcellweight,tlag,ncell)
        trans = trans0 + trans1

        trans2 = trans+np.transpose(trans)
        transP = []
        for t in trans2:
            if sum(t) != 0:
                transP.append(t/sum(t))
            else:
                transP.append(t)
        writeCSV(transP,of)
        eval,evec = np.linalg.eig(np.array(transP))
        print 'Implied Time Scale (unit iteration length): '+str(-tlag/np.log(eval[1]))


