#!/bin/bash
#SBATCH --job-name=364v2.3_MTHFR
#SBATCH --partition=dept_gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

# Load the NAMD module
module load namd/2.14/multicore-CUDA

# Set  variables
input_file=step5_production2
outputname=step5_2

# Restart the simulation
namd2 +idlepoll +p$SLURM_CPUS_PER_TASK +devices 0 ${input_file}.inp > ${outputname}.out

