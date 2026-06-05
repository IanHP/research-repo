#!/bin/bash

#SBATCH --job-name=cpptraj
#SBATCH --output=bash-op.txt

module purge
module load pixi

pixi run python create_export_cpp.py

echo "done"