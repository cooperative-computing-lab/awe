#!/usr/bin/env bash

home="$1"
gendir="$2"
run=$3
clone=$4
gen=$5

me=$(basename $0)


echo "[$me] $gendir"


complete='FALSE'

finalfile=$gendir/ala2.pdb

if [ -f $finalfile ]; then
    complete='TRUE'
    echo "[$me] $finalfile exists"
else
    complete='FALSE'
    echo "[$me] $finalfile does not exits. Listing directory contents:"
    ls $gendir
fi

if [ $complete == 'FALSE' ]; then
    exit 1
fi