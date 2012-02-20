import os
import sys
import numpy as np

def get_params(nw, ti):
    values = dict()
    values['REGEX_RUNS'] = nw
    values['REGEX_DATA'] = 'data.%s' % ti
    return values

def copy_replace(infilename, outfilename, values):
    """
    Copies file from infile to outfile line by line, replacing strings as it goes.
    Values is a dict mapping regexs to values.
    """
    infile = open(infilename)
    outfile = open(outfilename, "w")
    for ln in infile:
        for regex, value in values.iteritems():
            ln = ln.replace(regex, str(value))
        outfile.write(ln)
    infile.close()
    outfile.close()


if __name__ == "__main__":
    t0 =int( sys.argv[1]) # initial timestep
    ts =int( sys.argv[2]) # number of timesteps
#    nw =int( sys.argv[3]) # total number of walkers
    numstates = 100
    
#    tw = nw/numstates # walkers per state
    tw = 20 # we will increase from 3 to 5
    nw = int(np.loadtxt(os.path.join('startpdbs','busywalkers.txt')))

    for ti in range(t0,t0+ts):
#        assert not os.path.exists('data.%s' % ti)
        copy_replace('template-params.cfg', 'params.cfg', get_params(nw, ti))
        fawcmd = 'faw params.cfg'
        os.system(fawcmd)
        resamplecmd = 'python scripts/resample2.py %d' % (ti)
        os.system(resamplecmd)
        nw = int(np.loadtxt(os.path.join('startpdbs','busywalkers.txt')))

