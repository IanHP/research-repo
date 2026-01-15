#!/bin/bash
#SBATCH --job-name=WTv10_MTHFR
#SBATCH --partition=dept_gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

# Load the NAMD module
module load namd/2.14/multicore-CUDA

# Set initial variables
equi_prefix=step4_equilibration
prod_prefix=step5_production
prod_step=step5
cnt=1     # Initialize counter
cntmax=3    # Set the maximum count, adjust as needed

# Run equilibration
#namd2 +idlepoll +p10 ${equi_prefix}.inp > ${equi_prefix}.out
# namd2  ${equi_prefix}.inp > ${equi_prefix}.out
#namd2 +idlepoll +p$SLURM_TASKS_PER_NODE +devices 0 ${equi_prefix}.inp > ${equi_prefix}.out

namd2 +idlepoll +p$SLURM_CPUS_PER_TASK +devices 0 ${equi_prefix}.inp > ${equi_prefix}.out

# Set names for production run
inputname="${equi_prefix}"
outputname="${prod_step}_1"  # Assuming you want to keep this naming convention

# Prepare the input file for production run
# Change input and output names from template file
sed "s/${equi_prefix}/${inputname}/g" ${prod_prefix}.inp | \
    sed "s/${prod_prefix}/${outputname}/g" > ${prod_step}_run.inp

# Run the simulation for the production step
# namd2 +idlepoll +p$SLURM_TASKS_PER_NODE +devices 0 ${prod_step}_run.inp > ${outputname}.out

namd2 +idlepoll +p$SLURM_CPUS_PER_TASK +devices 0 ${prod_step}_run.inp > ${outputname}.out

# namd2 +idlepoll +p10 ${prod_step}_run.inp > ${outputname}.out
