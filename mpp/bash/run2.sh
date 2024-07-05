#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=8
#SBATCH --ntasks=256
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniele.corallo@kit.edu
#SBATCH --partition=cpuonly

alias mpirun="mpirun --bind-to core --map-by core --mca opal_warn_on_missing_libcuda 0"


mpirun M++ Experiment2 plevel=4 level=8 vtkplot=0 precision=10 logfile=log_exp2_uniform WeightedAssemble=0