import os
import shutil
import numpy as np

numstates = 100
nw = 20
colors = 2

def ndx(i,j,k):
    """k is color, i is state, j is walker """
    return i * nw + j + k * nw * numstates

if __name__ == '__main__':
    p=np.loadtxt('Populations.dat')
    dstdir='startpdbs'
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    else:
        cmd=('rm -fr startpdbs')        
        os.system(cmd)
    for k in range(colors):
        for i in range(numstates):
            for j in range(nw):
                idx = ndx(i,j,k)

		dst='weight%d.txt' % idx
		np.savetxt(os.path.join(dstdir,dst),[p[i] / nw] )
		dst='color%d.txt' % idx
		np.savetxt(os.path.join(dstdir,dst), [k])
		srcdir='/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs/'
		dst='state%d.pdb' % idx
		src='State%d-%d.pdb' % (i,j+k*nw)
		shutil.copy(os.path.join(srcdir,src),os.path.join(dstdir,dst))
		print 'src = %s, dst = %s' % (src,dst)
    np.savetxt(os.path.join(dstdir,'busywalkers.txt'),[colors*numstates*nw])
