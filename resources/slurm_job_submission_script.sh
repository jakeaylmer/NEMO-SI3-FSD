#!/bin/bash

# Example SLURM job submission script for NEMO-SI3 compiled with OpenMPI on the
# University of Reading Academic Computing Cluster 2 (UoR RACC2). This should
# be placed in a compiled configuration's run directory.

#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --nodes=1-1
#SBATCH --output=%x-%j.out
#SBATCH --time=1:0:0

module load mpi/openmpi-x86_64

# Do not use srun: the version of OpenMPI on RACC2 does not support it and
# mpiexec is recommended by OpenMPI anyway:
mpiexec ./nemo
