#!/bin/bash

#SBATCH --job-name=01category
#SBATCH --output=bash-op.txt
#SBATCH --time=0:05:00
#SBATCH --cluster=SMP

module purge
module load pixi

echo "start"

pixi run python -u categorize_interact.py

echo "done"
