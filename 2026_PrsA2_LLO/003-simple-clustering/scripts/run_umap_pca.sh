#!/bin/bash

#SBATCH --job-name=potn_inter
#SBATCH --output=bash-op.txt

module purge
module load pixi

pixi run python umap_graph.py
pixi run python pca_graph.py

echo "done"
