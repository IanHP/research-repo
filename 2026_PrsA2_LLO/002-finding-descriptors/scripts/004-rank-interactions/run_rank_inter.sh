#!/bin/bash

#SBATCH --job-name=004rank
#SBATCH --output=bash-op.txt
#SBATCH --time=4:00:00

module purge
module load pixi

echo "start"

pixi run python -u rank_interactions.py

echo "done"
