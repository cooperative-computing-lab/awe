""" resampling procedure for Adaptive Weighted Ensemble from Darve, E. (2011)
The objective is to maintain a number of walkers on each state and update the weights to accelerate convergence 
This is a 2-color algorithm"""

import numpy as np
import scipy.io
import itertools as it
import os
import shutil
import sys

trajname = 'discrete.traj'
outputname = 'ala2.pdb'
outputdir = 'output'
weightname = 'weight.txt'
colorname = 'color.txt'
datadir = 'Data'
gensname = 'Gens.lh5'
ndxname = 'AtomIndices.dat'
projname = 'ProjectInfo.h5'
startpdbs = 'startpdbs/'
fluxhistory = 'fluxhistory.txt'
transitionshistory = 'transitionshistory.txt'
FATALFILE = -1
BREAKPOINT = -2
#
WRITESTARTPDBS = False
# KEY VALUES
numstates = 100
workers = 20
sinkstates = [ [0,39,45,58,66,24,71], [6] ]
colors = 2
# END KEY VALUES
startstate = [] # list of start states for all walkers
endstate = [] # list of end states for all walkers
walkertopath = [] # list of end directory for all walkers
weights = [] # list of weights for all walkers
startcolor = []
endcolor = []
transitions = np.zeros((colors,colors))
fluxes=np.zeros((colors,colors))
fluxintoF=np.zeros(numstates)
fluxoutF=np.zeros(numstates)

def write_matrix():
	pass

def resample(weights,endstate,weighthistory,offset,mycolor):
        """ resample """
	print "resample color %s offset = %s" % (mycolor,offset)
        result=offset
        newweights=[]
        wr=np.arange(0,len(weights)) # walker range (cells x walkers / cell)
        er=np.arange(0,numstates) # state range (cells)
        eps = np.finfo(np.double).eps
        list1 = []
	try:
            if not os.path.exists(weighthistory):
           	file(weighthistory, 'w').close()    
	    wh = open(weighthistory, "a") 
	except IOError as e:
            print "IOError opening %s" % weighthistory

        for i in er:
            activewlk = 0
            ww = np.where(endstate==i)
            list0 = wr[ww]
            wi = weights[ww]
            ind = np.argsort(-wi)
            list0 = list(list0[ind])
            W = np.sum (wi)
            tw =  W / workers
#            print i, " = ", list0, " W = ", W, " tw = ", tw

            if len(list0)>0:
                result += 1
                x=list0.pop()
                while True:
                    Wx = weights[x]
                    if (Wx+eps >= tw):
                        r = int(np.floor( (Wx+eps) / tw ))
#                    print "x= ", x, "r= ", r, " Wx= ", Wx, " tw= ", tw
                        for item in it.repeat(x,r):
                            list1.append(item)
                            newweights.append(tw)
                        activewlk += r
#                    print "residual= ", Wx-r*tw
                        if activewlk < workers and Wx-r*tw+eps > 0.0:
#                        print "adding residual"
                            list0.append(x)
                            weights[x]=Wx-r*tw
                        if len(list0)>0:
                            x=list0.pop()
                        else:
                            break
                    else:
                        if len(list0)>0:
                            y = list0.pop()
                            Wy = weights[y]
                            Wxy = Wx + Wy
                            p=np.random.random()
                            if p < Wy / Wxy:
                                x = y
#                        print "chose ", x, " Wxy= ", Wxy
                            weights[x]=Wxy
	newweights=np.array(newweights)
	newweights /= np.sum(newweights)
        # new workers and new weights are in list1 and newweights
	if WRITESTARTPDBS:
		print "creating output files for timestep %d" % (timestep+1)

		assert(os.path.exists(startpdbs)) # need to make sure startpdbs exist before calling
		try:
			wh.write(str(timestep)+'\t ')
			for i in er:
				wh.write(str(np.sum(newweights[np.where(endstate[list1]==i)]))+'\t ')
			wh.write('\n')
			wh.close()
		except IOError as e:
			print "IOError writing to %s" % weighthistory
 
		for wlk, wt in enumerate(newweights):
			np.savetxt (os.path.join(startpdbs,'weight%d.txt' % (wlk+offset*workers)), [wt])
			np.savetxt (os.path.join(startpdbs,'color%d.txt' % (wlk+offset*workers)),[mycolor])
			srcpdb = os.path.join(walkertopath[list1[wlk]],outputname)
			dstpdb = os.path.join(startpdbs,'state%d.pdb' % (wlk+offset*workers))
			shutil.copy(srcpdb, dstpdb)
	return result

def newcolor(oldcolor,newstate):
    """ what should the new color of the worker be? """
    for i in range(colors):
        if newstate in sinkstates[i]: 
            transitions[oldcolor,i] += 1
            return i
    transitions[oldcolor,oldcolor] += 1
    return oldcolor

def process_discrete_traj(trajfile,weightfile,colorfile):
    """ estimate start and end states and colors """
    try:
        open(trajfile)
        traj = np.loadtxt(trajfile)
    except IOError as e:
        # this should not happen
        print "IOError opening %s " % trajfile
        sys.exit(FATALFILE)
    try:
        open(weightfile)
        w = np.loadtxt(weightfile)
        weights.append(w)
        startstate.append(int(traj[0])) #0-based
        endstate.append(int(traj[-1]))  #0-based
        walkertopath.append(top)
    except IOError as e:
        # this should not happen
        print "IOError opening %s " % weightfile
        sys.exit(FATALFILE)
    try:
        open(colorfile)
        color = int(np.loadtxt(colorfile))
        startcolor.append(color)
	newc=newcolor(color,traj[-1])
        endcolor.append(newc)
	fluxes[color,newc] += w
	folded = colors-1
	if (newc == folded) and (startstate[-1] not in sinkstates[colors-1]): # hack, assume folded state is always the last color
           print "flux into 6 from %s" % startstate[-1]
           fluxintoF[startstate[-1]] += w
	elif (color == folded) and endstate[-1] not in sinkstates[colors-1]:
           print "flux out of 6 to %s" % endstate[-1] 
           fluxoutF[endstate[-1]] += w
    except IOError as e:
        # this should not happen
        print "IOError opening %s " % colorfile
        sys.exit(FATALFILE)
        
    
def cluster_if_missing(top,datafile):
    """ cluster a trajectory if it failed to be clustered """
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    dst = os.path.join(datafile,gensname)
    shutil.copy(gensname,dst)
    dst = os.path.join(top,ndxname)
    shutil.copy(ndxname,dst)
    shutil.copy('state0.pdb',top)
    shutil.copy('with-env',top)
    shutil.copy('env.sh',top)
    msmdir = os.path.join(top,outputdir)
    projfile = os.path.join(top,projname)
    if os.path.exists(projfile):
        os.unlink(projfile)
    # change to current directory to execute msmb2
    os.chdir(top)
    if os.path.exists('Trajectories'):
        for root,dirs,files in os.walk('Trajectories'):
            for f in files: os.unlink(os.path.join(root, f))
            for d in dirs: shutil.rmtree(os.path.join(root, d))
        shutil.rmtree('Trajectories')
    withenv='./with-env'
    envsh='./env.sh'
    pdbfnm='state0.pdb'
    msmcmd='%s %s ConvertDataToHDF.py -s %s -I %s' % (withenv,envsh,pdbfnm,'output')
    os.system(msmcmd)
    msmcmd='%s %s Assign.py' % (withenv,envsh)
    os.system(msmcmd)
    msmcmd='%s %s ConvertAssignToText.py' % (withenv,envsh)
    os.system(msmcmd)
    os.chdir(currentdir)

if __name__ == "__main__":
    # First get the paths and the discrete trajectories to determine where workers ended up
    me = sys.argv[0]
    timestep = int(sys.argv[1])
#    workers  = int(sys.argv[2]) 
    startdir = 'data.%s' % (timestep)
    currentdir=os.getcwd()
    for top, dir, files in os.walk(startdir):
        for nm in files:
            if nm == outputname: 
                #final output exists -- does discrete trajectory?
                trajfile=os.path.join(top,datadir,trajname)
                weightfile=os.path.join(top,weightname)
                colorfile=os.path.join(top,colorname)
                if os.path.isfile(trajfile):
                    process_discrete_traj(trajfile,weightfile,colorfile)
                else:
                    # simulation ran but MSMBuilder didn't -- have to rerun it
                    print "Missing MSM for %s. Running it " % top
                    datafile=os.path.join(top,datadir)
                    cluster_if_missing(top,datafile)
                    process_discrete_traj(trajfile,weightfile,colorfile)
    try:
	    if not os.path.exists(fluxhistory):
		    file(fluxhistory, 'w').close()    
	    fh = open(fluxhistory, "a")
    except IOError as e:
	    print "IOError opening %s" % fluxhistory
    try:
	    fh.write(str(timestep)+'\t ')
	    for i in range(colors):
		    for j in range(colors):
			    fh.write(str(fluxes[i,j])+'\t')
	    fh.write('\n')
	    fh.close()
    except IOError as e:
	    print "IOError writing %s" % fluxhistory
    
    try:
	    if not os.path.exists(transitionshistory):
		    file(transitionshistory, 'w').close()    
	    th = open(transitionshistory, "a")
    except IOError as e:
	    print "IOError opening %s" % transitionshistory
    try:
	    th.write(str(timestep)+'\t ')
	    for i in range(colors):
		    for j in range(colors):
			    th.write(str(transitions[i,j])+'\t')
	    th.write('\n')
	    th.close()
    except IOError as e:
	    print "IOError writing %s" % transitionshistory
    try:
	    fh = open('fluxintoF.txt','w')	    
	    np.savetxt('fluxintoF.txt',fluxintoF)
	    fh = open('fluxoutF.txt','w')	 
	    np.savetxt('fluxoutF.txt',fluxoutF)   
    except IOError as e:
	    print "IOError writing %s" % 'fluxintoF.xtx'

    endcolor = np.array(endcolor)
    if WRITESTARTPDBS:
	    if os.path.exists(startpdbs):
		    cmd=('rm -fr startpdbs')        
		    os.system(cmd)
	    os.makedirs(startpdbs)
    
    weights=np.array(weights)
    endstate=np.array(endstate)
    offset=0
    for cc in range(colors):
        indxs = np.where(endcolor == cc)
        colorweights=weights[indxs]
        colorendstate=endstate[indxs]
        weighthistory = 'weighthistory.%s.txt' % cc
        offset=resample(colorweights,colorendstate,weighthistory,offset,cc)
	print "after resample, offset = %s" % offset
    if WRITESTARTPDBS:
	    np.savetxt (os.path.join(startpdbs,'busywalkers.txt'), [offset*workers])  


            


    
                
                


        

    
    
                
    
            
            
            
    

    


