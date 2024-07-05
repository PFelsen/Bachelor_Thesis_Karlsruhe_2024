#!/bin/bash
mfile='-x LD_LIBRARY_PATH '
if [ $# -ge 0 ]
then
  mfile+='--machinefile '${1}
  echo "arg: $mfile"
fi

mpirun $mfile -n 128 M++ Experiment1 plevel=6 level=8 vtkplot=0 logfile=log/log_exp1_uniform precision=10
mpirun $mfile -n 128 M++ Experiment1 run=adaptive plevel=6 level=6 vtkplot=0 theta=0.2 print_eta_dist=1 logfile=log/log_exp1_adaptive6 precision=10 refinement_steps=2

mpirun $mfile -n 128 M++ Experiment2 plevel=6 level=8 vtkplot=0 logfile=log/log_exp2_uniform precision=10
mpirun $mfile -n 128 M++ Experiment2 run=adaptive plevel=6 level=6 vtkplot=0 theta=0.2 print_eta_dist=1 logfile=log/log_exp2_adaptive6 precision=10 refinement_steps=2

mpirun $mfile -n 128 M++ Experiment3 plevel=6 level=7 vtkplot=0 logfile=log/log_exp3_uniform precision=10
mpirun $mfile -n 128 M++ Experiment3 run=adaptive plevel=6 level=6 vtkplot=0 theta=0.2 print_eta_dist=1 logfile=log/log_exp3_adaptive6 precision=10 refinement_steps=2