#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --nodes=16
#SBATCH --ntasks=1024
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniele.corallo@kit.edu
#SBATCH --partition=cpuonly

COMMON_PARAM="vtkplot=0 precision=10"

alias mpirun="mpirun --bind-to core --map-by core --mca opal_warn_on_missing_libcuda 0"

mpirun M++ Experiment1 $COMMON_PARAM plevel=6 level=8 logfile=log/log_exp1_uniform
mpirun M++ Experiment1 $COMMON_PARAM run=adaptive plevel=6 level=6 theta=0.2 print_eta_dist=1 logfile=log/log_exp1_adaptive6 refinement_steps=2

mpirun M++ Experiment2 $COMMON_PARAM plevel=6 level=8 logfile=log/log_exp2_uniform
mpirun M++ Experiment2 $COMMON_PARAM run=adaptive plevel=6 level=6 theta=0.2 print_eta_dist=1 logfile=log/log_exp2_adaptive6 refinement_steps=2

mpirun M++ Experiment3 $COMMON_PARAM plevel=6 level=7 logfile=log/log_exp3_uniform
mpirun M++ Experiment3 $COMMON_PARAM run=adaptive plevel=6 level=6 theta=0.2 print_eta_dist=1 logfile=log/log_exp3_adaptive6 refinement_steps=2



