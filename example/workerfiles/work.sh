#!/usr/bin/env bash


echo "ls"
ls
echo

echo "creating result files"
cp -v structure.pdb structure2.pdb
echo '0' > cell2.dat
echo

echo "compressing to results.tar"
tar cvf results.tar structure2.pdb *.dat
echo


exit 0
