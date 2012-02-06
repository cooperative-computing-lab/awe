#!/usr/bin/env bash


gendir="$1"
run=$2
clone=$3
gen=$4

me=$(basename $0)

echo "[$me] $gendir $run $clone $gen"
cd $gendir

echo "[$me] Current files:"
ls

echo "[$me] extracting results:"
tar xvf result.tar.bz2 || exit 250

echo "[$me] removing tarfile:"
rm -v result.tar.bz2
