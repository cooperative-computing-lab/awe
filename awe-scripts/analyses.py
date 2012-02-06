""" Analyses of data produced by the Adaptive Weighted Ensemble
"""
import numpy as np
import scipy.io
import os
import sys

#import pylab
#import matplotlib 
#import matplotlib.pyplot

def populations(weightsfile='weighthistory_mod.txt',plotthem=True,WINDOW=5,NUMSTATES=100,IGNORE=0,listofstates=[]):
    results = []
    weightings=np.repeat(1.0, WINDOW) / WINDOW	
    wts=np.loadtxt(weightsfile)
    hi=wts.transpose()
    means = np.empty(NUMSTATES)
    devs  = np.empty(NUMSTATES)
    for i in range(1,NUMSTATES+1):
        hi1=np.convolve(hi[i,:], weightings)[WINDOW-1:-(WINDOW-1)]
        means[i-1]=np.mean(hi1[IGNORE:])
        devs[i-1]=np.std(hi1[IGNORE:])
        if plotthem and (i-1) in listofstates:
            results.append(hi1)
#            matplotlib.figure.Figure()
#            matplotlib.pyplot.plot(hi1,label='State %d' % i)
#    matplotlib.backends.backend_wx.Show()
    return means, devs, results

def ratematrix(T):
    """ compute populations from rate matrix """
    T=T.tocsr()
    n=T.shape[0] #number of states
    r=scipy.matrix(np.zeros((n,n)))
    r[0,:]=np.ones(n) # probability normalization
    b=np.zeros(n)
    b[0]=1.
    for i in range(1,n):
        for j in range(n):
            if i!=j:
                r[i,j] += T[j,i] # flux in
                r[i,i] -= T[i,j] # flux out
    p=np.linalg.linalg.solve(r,b)
    return p

def computeflow_i(numstates, dirnum):
    """ compute flux into cells
    using sum_i (w_i * N_ij / Ti) where the i are walkers with weight w_i 
    that crossed into state j N_ij times 
    with lifetime Ti before crossing into j 
    """
    dir = 'data.%d' % dirnum

    matrices = []
    weights = np.zeros((numstates+1,1))
    N = np.zeros((numstates+1,numstates+1))
    trajname = 'discrete.traj'
    FnTUnSym = "tCounts.UnSym.mtx"
    # get all discrete trajectories and transition count matrices for all RUNs except j
    for top, dirs, files in os.walk(dir):
        for nm in files:
            if nm == (trajname):
                traj=np.loadtxt(os.path.join(top,nm))
                w=np.loadtxt(os.path.join(top,"../weight.txt"))
                weights[traj[0]]=w
                t=scipy.io.mmread(os.path.join(top,FnTUnSym))
                matrices.append(t)
    # compress matrices
    n=np.zeros((numstates+1,numstates+1))
    for i in range(len(matrices)):
        N += matrices[i]
    ww = weights.transpose()
    NT = N.transpose()
    flowin = np.dot(ww,N) - ww*np.diag(N)
    flowout = np.dot(ww,NT) - ww*np.diag(NT)
    return flowin, flowout


    

    


        
                     
    
