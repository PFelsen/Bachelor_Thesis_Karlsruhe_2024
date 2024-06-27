# Bachelor_Thesis_Karlsruhe_2024
This repository contains the code for my Bachelor's thesis, completed in July 2024 at the Karlsruhe Institute of Technology. The Python folder contains my implementation of sparse grid quadrature, a numerical quadrature algorithm based on hierarchical interpolation for approximating integrals over 2- and 3-dimensional rectangles. The mpp folder contains a copy of the M++ software (https://www.math.kit.edu/ianm3/seite/mplusplus/en).

M++ (Meshes, Multigrid, and more) is a parallel finite element software developed at the Institute for Applied and Numerical Mathematics in the Scientific Computing group. It is used in various research projects and for teaching purposes. M++ is implemented in the C++ programming language and utilizes Open MPI to distribute the computational load across several computing nodes. The library TASMANIAN (https://github.com/ORNL/TASMANIAN) is used for stochastic collocation.

To reproduce the numerical results in my Bachelor's thesis, refer to the Jupyter Notebook located at mpp/notebooks/SC/UQ_SC_MC_Laplace.ipynb. For execution instructions, please refer to the M++ README file.
