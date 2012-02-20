#!/afs/crc.nd.edu/user/c/cabdulwa/Public/apps/python/x86_64/2.7.1/bin/python
import numpy
import scipy.io
from msmbuilder import Serializer,MSMLib

if __name__ == "__main__":
    print """ Transforming assignments into a discrete trajectory
          and count matrix """

    FnTUnSym ="Data/tCounts.UnSym.mtx"
    NumStates=101 # one more than there are states
    Assignments = Serializer.LoadData("Data/Assignments.h5")

    Counts=MSMLib.GetCountMatrixFromAssignments(Assignments, NumStates, LagTime=1, Slide=True)
    scipy.io.mmwrite(FnTUnSym, Counts) # ignore the NumStates, NumStates transition

    Assignments = Assignments.flatten()
    numpy.savetxt ("Data/discrete.traj", Assignments, fmt="%-1d")
