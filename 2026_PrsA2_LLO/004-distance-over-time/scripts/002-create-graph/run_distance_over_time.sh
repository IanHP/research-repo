#!/bin/bash

#SBATCH --job-name=004rank
#SBATCH --output=bash-op.txt
#SBATCH --time=0:05:00
#SBATCH --cluster=SMP
#SBATCH --partition=smp



module purge
module load pixi

echo "start"

pixi run python -u distance_over_time.py

echo "done"
