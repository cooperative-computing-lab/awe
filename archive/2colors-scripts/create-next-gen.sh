#!/usr/bin/env bash

set -o errexit

oldgendir="$1"
newgendir="$2"
run=$3
clone=$4
gen=$5
newgen=$6

me=$(basename $0)

real_newgendir=$(readlink -f $newgendir)

echo "[$me] $oldgendir $newgendir $run $clone $gen $newgen"

[ ! -d $newgendir ] &&
mkdir -pv $newgendir

echo "[$me] copying files from G$gen"
cp -v $oldgendir/{*.conf,*.inp,*.psf} $newgendir
cp -v $oldgendir/final.pdb $newgendir/start.pdb

echo "[$me] creating new tprfile"
cd $oldgendir

