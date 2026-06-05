#!/bin/bash

#SBATCH --job-name=potn_inter
#SBATCH --output=bash-op.txt

module purge
module load pixi

pixi run python all_potential_interactions.py

echo "done"
