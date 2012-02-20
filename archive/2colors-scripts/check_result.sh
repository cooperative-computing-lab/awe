#!/usr/bin/env bash


gendir=$1
run=$2
clone=$3
gen=$4

me=$(basename $0)

logfile=$gendir/my_md_logfile.txt

echo "[$me] $gendir $run $clone $gen"


#echo "[$me] checking if simulations aborted:"
#result=$(grep ABORT $logfile | sed "s: *::")
#echo "    found: |$result|"

#if [ ! "$result" == "" ]; then
#        exit 1
#fi


#echo "[$me] Checking for Segfault"
#result=$(grep -i 'segmentation fault' $logfile | sed "s: *::")
#echo "    found: |$result|"

#if [ ! "$result" == "" ]; then
#        exit 1
#fi


echo "[$me] Checking for SUCCESS"
result=$(tail $logfile | grep SUCCESS | sed "s: *::")
echo "    found: |$result|"

if [ "$result" == "" ]; then
        exit 1
fi
