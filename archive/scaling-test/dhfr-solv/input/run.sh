#!/usr/bin/env bash

set -o errexit

echo "directory contents"
ls
echo

echo "setup and load modules"
source ~cabdulwa/.bash_modules
module load gromacs/4.5.3

### disable backups
export GMX_MAXBACKUP=-1
echo

echo "running simulation"
mdrun -s *.tpr -nt 1
echo

echo "setup result files"
cp -v structure.pdb structure2.pdb
cp -v weight.dat weight2.dat
cp -v color.dat color2.dat
cp -v cell.dat cell2.dat
echo

echo "compressing"
tar cvf results.tar structure2.pdb *.dat
echo

exit 0