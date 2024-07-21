# Bachelor_Thesis_Karlsruhe_2024
This repository contains the code for my Bachelor's thesis, completed in July 2024 at the Karlsruhe Institute of Technology. The Python folder contains my implementation of sparse grid quadrature, a numerical quadrature algorithm based on hierarchical interpolation for approximating integrals over 2- and 3-dimensional rectangles. 

To reproduce the numerical results in my Bachelor's thesis, refer to the Jupyter Notebook located at mpp/notebooks/SC/UQ_SC_MC_Laplace.ipynb. For execution instructions, please refer to the M++ README file.


## Repository Structure

### Folders

- **mpp**: The mpp folder contains a copy of the M++ software (https://www.math.kit.edu/ianm3/seite/mplusplus/en). This version will not be changed again. It is a copy of the branch 767-repair-stochasticcollocation, and there is no version control over it through the main M++ Git repository (https://gitlab.kit.edu/kit/mpp/mpp).
For execution instructions, please refer to the M++ README file.
M++ (Meshes, Multigrid, and more) is a parallel finite element software developed at the Institute for Applied and Numerical Mathematics in the Scientific Computing group. It is used in various research projects and for teaching purposes. M++ is implemented in the C++ programming language and utilizes Open MPI to distribute the computational load across several computing nodes. The library TASMANIAN (https://github.com/ORNL/TASMANIAN) is used for stochastic collocation.
Keep in mind that the current code can only be run on a single core. The parallel implementation has not yet been implemented due to the negative weights in the quadrature rule used.

- **MonteCarlo_logfiles**: Logfiles for Monte Carlo calculations with 100,000 samples.
- **python**: Python scripts for testing and comparing integration methods.

### Key Files

- **mpp/notebooks/SC/UQ_SC_MC_Laplace.ipynb**: Jupyter Notebook for reproducing thesis results.
- **python/sparse_grid_integration/sparse_grid_integration_2D.py**: Implementation of sparse grid integration for 2D functions.
- **python/sparse_grid_integration/sparse_grid_integration_3D.py**: Implementation of sparse grid integration for 3D functions.
- **python/sparse_grid_integration/sparse_grid_vs_full_grid_quadrature_2D.py**: Comparison of sparse grid and full grid quadrature in 2D.
- **python/sparse_grid_integration/sparse_vs_full_error_polynomial_plots.py**: Error plots of sparse vs. full grid quadrature.
- **python/sparse_grid_integration/sparse_vs_full_grid_random_polynomial_test.py**: Tests comparing sparse and full grid methods with random polynomial functions.
- **python/convergence_tests/convergence_tests_kernel.py**: Convergence tests for different kernels.
- **python/convergence_tests/convergence_tests_standard.py**: Convergence tests for standard kernel.
- **python/convergence_tests/DataFrames**: Contains the raw DataFrames created with the Jupyter Notebook to calculate the convergence results.

## Requirements

- Python 3.x
- Jupyter Notebook
- M++ software (included in the repository)
