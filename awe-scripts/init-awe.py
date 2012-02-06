import os
import shutil
import numpy as np

if __name__ == '__main__':
	p=np.loadtxt('Populations.dat')
	numstates = 100
	nw = 10
	for i in range(numstates):
		for j in range(nw):
			dstdir='./startpdbs'
			dst='weight%d.txt' % (i * nw + j)
			np.savetxt(os.path.join(dstdir,dst),[p[i]])
                        srcdir='/afs/crc.nd.edu/user/i/izaguirr/Public/ala2/faw-protomol/PDBs/'
                        dst='state%d.pdb' % (i * nw + j)
                        src='State%d-%d.pdb' % (i,j)
			shutil.copy(os.path.join(srcdir,src),os.path.join(dstdir,dst))

