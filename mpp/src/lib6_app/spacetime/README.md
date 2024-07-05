To build spacetime with {1,2,3}-D space-cells use

- `cmake .. -DUSE_SPACETIME=ON -DSPACE_DIM=1`
- `cmake .. -DUSE_SPACETIME=ON -DSPACE_DIM=2`
- `cmake .. -DUSE_SPACETIME=ON -DSPACE_DIM=3`

For Horeka-Users the following modules have to be loaded:

-`module load compiler/gnu/8 mpi/openmpi/4.1 devel/cmake/3.23`

For the PDE-Cluster you need:

-`module load OpenMPI/4.1.4-GCC-11.3.0 CMake/3.23.1-GCCcore-11.3.0 ImageMagick/7.1.0-37-GCCcore-11.3.0 Valgrind/3.18.1-gompi-2022a`