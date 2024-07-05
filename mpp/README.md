## What is M++?
A C++ library providing Finite Element Methods (FEM) to solve Partial Differential Equations (PDE).
It uses MPI to distribute the computations on parallel processes and nodes.

The gitlab-group https://gitlab.kit.edu/kit/mpp provides an overview of associated projects. 


### Prerequisites (build environment)
M++ uses:
* CMake (https://cmake.org/download/)  
* GCC >= 11 (C++ 20 needs to be fully supported)   
* BLAS (Basic Linear Algebra Subroutines https://www.netlib.org/blas/)  
* LAPACK (Linear Algebra PACKage https://www.netlib.org/lapack/)  
* Open MPI >= 4.0 (https://www.open-mpi.org/)  

Optional:  
* FFTW   

The following library are automatically provided via submodules:
* TASMANIAN  
* Google Benchmark  
* GoogleTest  
* pugixml  
* SPRNG5  
* SuperLU  

### Installing M++
* You can install the packages required to build M++ on Ubuntu 22.04 with the following command:<br>
```sudo apt-get install git cmake build-essential libblas-dev liblapack-dev libopenmpi-dev libtirpc-dev gcc```  
  
* To clone the project with the necessary submodules run:  
```git clone --recurse-submodules https://gitlab.kit.edu/kit/mpp/mpp```

* Alternatively, clone it and update the submodules by hand:  
```git clone https://gitlab.kit.edu/kit/mpp/mpp```  
```cd mpp```  
```git submodule update --init --recursive```   

* To build M++, run the following commands in the ```mpp``` root directory, where <x> represents the number of processes to compile with:  
```mkdir build```  
```cd build```     
```cmake ..```  
```make -j <x>```

* Alternatively, simply run:                 
```./init.sh```

### Run tests
* If built with ```BUILD_TESTS=ON```, Unit tests can be executed via (provided a python3-environment with installed ```matplotlib```):                           
```cd build```                                                                  
```ctest```                                                                     
```python3 mppyrun.py --mpi_tests=1 --mute=1```                                         

### Important build arguments   
| build argument | parameters | description | 
| :--- | :--- | :--- |
| ```BUILD_TESTS``` | ```ON``` , ```OFF``` | Builds the core test-library  |
| ```BUILD_TUTORIAL``` | ```ON``` , ```OFF``` | Enables the use of the included tutorial-module (e.g. acoustic and cransport siscretizations and Assembles)  |
| ```BUILD_TUTORIAL_TESTS```| ```ON``` , ```OFF``` |  Builds the tutorial-test-library  |
| ```USE_SPACETIME``` | ```ON``` , ```OFF``` | Enables the use of the included spacetime-module (also sets TIME_DIM=1, hence each Point has a time-component)  |
| ```BUILD_UQ```| ```ON``` , ```OFF``` | Builds the included uncertainty quantification-module (e.g. (Budgeted) MultiLevelMonteCarlo, SGD for StochasticOptimization)  |
| ```SPACE_DIM``` | ```1``` , ```2```,  ```3``` | Sets the number of dimensions for a Point-object. (needs to be at least the dimension required for given problem) | 
| ```MPP_BUILD_TYPE```| ```MppRelease```, ```MppDebugFast```, ```MppDebug``` | Enables different types of compiler-optimization. For Debuggig use ```MppDebug```|  


Note: These build arguments are used by calling cmake in the following way, for example in the build-directory:  
```cmake .. -DUSE_SPACETIME=ON -DBUILD_TUTORIAL=ON```

If you used another build-configuration, please delete the automatically created `CMakeCache.txt`.

