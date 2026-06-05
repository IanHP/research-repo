#!/bin/bash

#SBATCH --job-name=cpptraj
#SBATCH --output=bash-op.txt

module purge
module load amber/24

cpptraj ph5_update.cpp
echo "done"
