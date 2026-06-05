#!/bin/bash

#SBATCH --job-name=cpptraj
#SBATCH --output=bash-op.txt

module purge
module load amber/24

cpptraj cpp/export_prsa2-llo-ph5-r1_op.cpp
cpptraj cpp/export_prsa2-llo-ph5-r2_op.cpp
cpptraj cpp/export_prsa2-llo-ph5-r3_op.cpp

cpptraj cpp/export_prsa2-llo-ph7-r1_op.cpp
cpptraj cpp/export_prsa2-llo-ph7-r2_op.cpp
cpptraj cpp/export_prsa2-llo-ph7-r3_op.cpp
echo "done"
