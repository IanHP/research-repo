#!/bin/bash

#SBATCH --job-name=distance
#SBATCH --output=bash-op.txt
#SBATCH --time=4:00:00

module purge
module load pixi

echo "start"

pixi run python -u each_frame_interaction_distance.py

echo "done"
