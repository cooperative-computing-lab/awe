#!/usr/bin/env bash

source ~cabdulwa/.bash_modules
module use ~cabdulwa/Public/modulefiles
module load gromacs/4.0.7

pdb=$1

python prepare_4.0.7.py $pdb input/RUN10.nocaps.top input/run.mdp
