""" resampling procedure for Adaptive Weighted Ensemble from Darve, E. (2011)
The objective is to maintain a number of walkers on each state and update the weights to accelerate convergence """

import numpy as np
import scipy.io
import itertools as it
import os
import shutil
import sys

def cluster_if_missing():
    """ cluster a trajectory if it failed to be clustered """
    

if __name__ == "__main__":
    # First get the paths and the discrete trajectories to determine where workers ended up
    trajname = 'discrete.traj'
    outputname = 'ala2.pdb'
    outputdir = 'output'
    weightname = 'weight.txt'
    weighthistory = 'weighthistory.txt'
    datadir = 'Data'
    gensname = 'Gens.lh5'
    ndxname = 'AtomIndices.dat'
    projname = 'ProjectInfo.h5'
    FATALFILE = -1
    BREAKPOINT = -2
    missingmsm = False
    numstates = 100
    workers = 20
    me = sys.argv[0]
    timestep = int(sys.argv[1])
#    workers  = int(sys.argv[2]) 

    startdir = 'data.%s' % timestep
    startstate = [] # list of start states for all walkers
    endstate = [] # list of end states for all walkers
    walkertopath = [] # list of end directory for all walkers
    weights = [] # list of weights for all walkers
    currentdir=os.getcwd()
    for top, dir, files in os.walk(startdir):
        for nm in files:
            if nm == outputname: 
                #final output exists -- does discrete trajectory?
                trajfile=os.path.join(top,datadir,trajname)
                weightfile=os.path.join(top,weightname)
                if os.path.isfile(trajfile):
                    try:
                        open(trajfile)
                        traj = np.loadtxt(trajfile)
                    except IOError as e:
                        # this should not happen
                        print "Fatal file error opening %s " % trajfile
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
                        print "Fatal file error opening %s " % weightfile
                        sys.exit(FATALFILE)
                else:
                    # simulation ran but MSMBuilder didn't -- have to rerun it
                    print "Missing MSM for %s. Running it " % top
                    datafile=os.path.join(top,datadir)
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
#                    sys.exit(BREAKPOINT)
    if not missingmsm:
        weights=np.array(weights)
        newweights=[]
        endstate=np.array(endstate)
        wr=np.arange(0,len(weights)) # walker range (cells x walkers / cell)
        er=np.arange(0,numstates) # state range (cells)
        eps = np.finfo(np.double).eps
        list1 = []
        busycells = 0

        if not os.path.exists(weighthistory):
            file(weighthistory, 'w').close()    
        wh = open(weighthistory, "a") 

        for i in er:
            activewlk = 0
            ww = np.where(endstate==i)
            list0 = wr[ww]
            wi = weights[ww]
            ind = np.argsort(-wi)
            list0 = list(list0[ind])
            W = np.sum (wi)
            tw =  W / workers
            print i, " = ", list0, " W = ", W, " tw = ", tw

            if len(list0)>0:
                busycells += 1
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
        # normalize newweights -- roundoff can make their sum different than 1
        newweights=np.array(newweights)
        newweights /= np.sum(newweights)
        # new workers and new weights are in list1 and newweights
        print "creating output files for timestep %d" % (timestep+1)
        startpdbs = 'startpdbs/'
        if not os.path.exists(startpdbs):
            file(startpdbs,'w').close()
        for root,dirs,files in os.walk(startpdbs):
            for f in files: os.unlink(os.path.join(root, f))
            for d in dirs: shutil.rmtree(os.path.join(root, d))
        # need to write weight%walker.txt with newweights
        # need to write state%walker.xyz with positions
        wh.write(str(timestep)+'\t ')
        for i in er:
            wh.write(str(np.sum(newweights[np.where(endstate[list1]==i)]))+'\t ')
        wh.write('\n')
        wh.close()
        np.savetxt (os.path.join(startpdbs,'busywalkers.txt'), [busycells*workers])  
        for wlk, wt in enumerate(newweights):
            np.savetxt (os.path.join(startpdbs,'weight%d.txt' % wlk), [wt])
            srcpdb = os.path.join(walkertopath[list1[wlk]],outputname)
            dstpdb = os.path.join(startpdbs,'state%d.pdb' % wlk)
            shutil.copy(srcpdb, dstpdb)
        



            


    
                
                


        

    
    
                
    
            
            
            
    

    


