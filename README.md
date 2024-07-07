# Bachelor_Thesis_Karlsruhe_2024
This repository contains the code for my Bachelor's thesis, completed in July 2024 at the Karlsruhe Institute of Technology. The Python folder contains my implementation of sparse grid quadrature, a numerical quadrature algorithm based on hierarchical interpolation for approximating integrals over 2- and 3-dimensional rectangles. 

To reproduce the numerical results in my Bachelor's thesis, refer to the Jupyter Notebook located at mpp/notebooks/SC/UQ_SC_MC_Laplace.ipynb. For execution instructions, please refer to the M++ README file.


# Folders
## mpp
The mpp folder contains a copy of the M++ software (https://www.math.kit.edu/ianm3/seite/mplusplus/en). This version will not be changed again. It is a copy of the branch 767-repair-stochasticcollocation, and there is no version control over it through the main M++ Git repository.
For execution instructions, please refer to the M++ README file.
M++ (Meshes, Multigrid, and more) is a parallel finite element software developed at the Institute for Applied and Numerical Mathematics in the Scientific Computing group. It is used in various research projects and for teaching purposes. M++ is implemented in the C++ programming language and utilizes Open MPI to distribute the computational load across several computing nodes. The library TASMANIAN (https://github.com/ORNL/TASMANIAN) is used for stochastic collocation.
Keep in mind that the current code can only be run on a single core. The parallel implementation has not yet been implemented due to the negative weights in the quadrature rule used.

## MonteCarlo_logfiles
Contains the logfiles for Monte Carlo (MC) calculations using 100,000 samples. These are used as reference values to validate the Stochastic Collocation method.

## python

This folder contains Python programs to test and compare the integration of various mathematical functions using sparse grid and full grid quadrature methods. The scripts provided are designed to handle integration in both two-dimensional and three-dimensional spaces.

### File Descriptions

### `sparse_grid_integration_2D.py`
This file contains the implementation of sparse grid integration for two-dimensional functions. It includes functions for defining piecewise linear basis functions, computing hierarchical surpluses, and performing the sparse grid integration.

#### Main Functions:
- `phi(x)`: Computes the value of the piecewise linear function.
- `phi_j_k(x, j, k, a, b)`: Computes the value of the scaled and translated piecewise linear function.
- `int_phi(l, a, b)`: Computes the integral values for the function.
- `delta(f, l_1, l_2, i_1, i_2, a, b, c, d)`: Computes the hierarchical surplus for a given function at specific levels and indices.
- `J(k)`: Generates the index set for a given level.
- `sparse_grid_integration(f, L, X, Y)`: Computes the integral for a given function in two dimensions over \([a, b] \times [c, d]\).
- `function_test()`: Tests the sparse grid integration with random polynomial, exponential, and trigonometric functions in two dimensions.

### `sparse_grid_integration_3D.py`
This file contains the implementation of sparse grid integration for three-dimensional functions. It includes functions for defining piecewise linear basis functions, computing hierarchical surpluses, and performing the sparse grid integration in three dimensions.

#### Main Functions:
- `phi(x)`: Computes the value of the piecewise linear function.
- `phi_j_k(x, j, k, a, b)`: Computes the value of the scaled and translated piecewise linear function.
- `int_phi(l, a, b)`: Computes the integral values for the function.
- `delta(f, l_1, l_2, l_3, i_1, i_2, i_3, a, b, c, d, e, g)`: Computes the hierarchical surplus for a given function at specific levels and indices in three dimensions.
- `J(k)`: Generates the index set for a given level.
- `sparse_grid_integration(f, L, X, Y, Z)`: Computes the integral for a given function in three dimensions over \([a, b] \times [c, d] \times [e, g]\).
- `function_test()`: Tests the sparse grid integration with random polynomial, exponential, and trigonometric functions in three dimensions.

### `sparse_grid_vs_full_grid_quadrature_2D.py`
This file contains scripts to compare the performance and error of sparse grid quadrature versus full grid quadrature in two dimensions. It includes functions to compute the integral using both sparse and full grids and plots the errors.

### `sparse_vs_full_error_polynomial_plots.py`
This file includes scripts to plot the errors of sparse grid and full grid quadrature methods for polynomial functions of varying degrees. It showcases the efficiency and accuracy of sparse grid methods in higher dimensions.

### `sparse_vs_full_grid_random_polynomial_test.py`
This file contains tests comparing the performance of sparse grid and full grid methods using randomly generated polynomial functions. It evaluates the accuracy and computational time of both methods across different polynomial degrees and levels of refinement.
