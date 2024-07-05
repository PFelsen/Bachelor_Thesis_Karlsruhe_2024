#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=14
#SBATCH --ntasks=1024
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniele.corallo@kit.edu
#SBATCH --partition=cpuonly

alias mpirun="mpirun --bind-to core --map-by core --mca opal_warn_on_missing_libcuda 0"

mpirun M++ riemannDoubleSquare normal_y=0.0 run=single time_deg=2 space_deg=2 MeshesVerbose=2 MeshVerbose=1 ConfigVerbose=4 ExcludedResults="L2SpaceAtT, L2SpaceAtT_int, ||u_proj-u_h||_DG, DiscNorm, EE, Conf" plevel=4 level=7 logfile=log/scaling${SLURM_NTASKS}.txt

