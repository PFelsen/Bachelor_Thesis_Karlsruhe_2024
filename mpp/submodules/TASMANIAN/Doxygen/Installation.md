# Installation

[TOC]

### Document Sections
* Requirements
* Install using CMake: the preferred way
* Install with (basic) GNU Make
* Install with Python Pip
* Install with Spack
* Install on MS Windows platform
* Install folder structure
* Linking to Tasmanian: CMake Package Config
* Known issues

### Requirements

Minimum requirements to use Tasmanian:
* a C/C++ compiler and either [CMake](https://cmake.org/) or [GNU Make](https://www.gnu.org/software/make/) build engines.

Recommended additional features:
* [Python](https://www.python.org/) with with [NumPy](http://www.numpy.org/) 1.10 (or newer) and [CTypes](https://docs.python.org/3/library/ctypes.html) packages
* [Basic Linear Algebra Subroutine (BLAS)](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) implementation
* [Linear Algebra PACKage](http://www.netlib.org/lapack/) implementation
* [OpenMP](https://en.wikipedia.org/wiki/OpenMP) implementation (usually included with the compiler)

Optional features:
* Acceleration using Nvidia [linear algebra libraries](https://developer.nvidia.com/cublas) and custom [CUDA kernels](https://developer.nvidia.com/cuda-zone)
* Acceleration using AMD ROCm [linear algebra libraries](https://rocsparse.readthedocs.io/en/master/) and custom [HIP kernels](https://rocmdocs.amd.com/en/latest/ROCm_API_References/HIP-API.html)
* Acceleration using Intel OneAPI [oneMKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html) and custom [DPC++ kernels](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html)
* GPU out-of-core algorithms using the [UTK MAGMA library](http://icl.cs.utk.edu/magma/)
* Basic [Python matplotlib](https://matplotlib.org/) support
* Fully featured [MATLAB/Octave](https://www.gnu.org/software/octave/) interface via wrappers around the command-line tool
* Fortran 2003 interface using [gfortran](https://gcc.gnu.org/wiki/GFortran), [flang](https://flang.llvm.org/docs/), [ifort](https://software.intel.com/en-us/intel-compilers), [pgfortran](https://www.pgroup.com/index.htm), [xlf2003](https://www.ibm.com/products/xl-fortran-linux-compiler-power)
* Addon templates for the [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface)
* [Doxygen](http://www.doxygen.org/) documentation

**Note:** with the exception of the MAGMA library and some of the more exotic compilers, the rest of the software is included in most Linux distributions (or HomeBrew) or there are vendor supported software repositories.

| Feature | Tested versions     | Recommended      |
|----|----|----|
| gcc     | 9 - 11              | any              |
| clang   | 10 - 15             | any              |
| icc     | 18.1                | 18.1             |
| xl      | 16.1                | 16.1             |
| pgi     | 19.10 - 20.4        | 20.4             |
| cmake   | 3.19 - 3.25         | 3.22             |
| python  | 3.8 - 3.10          | any              |
| anaconda| 5.3                 | 5.3              |
| OpenBlas| 0.2.18 - 3.08       | any              |
| ESSL    | 6.2                 | 6.2              |
| CUDA    | 10.0 - 12           | 12.1             |
| ROCm    | 5.5 -5.7            | any              |
| DPC++   | 2021 - 2023         | 2023             |
| libiomp | 5.0                 | 5.0              |
| MAGMA   | 2.5.1 - 2.6.1       | 2.6.1            |
| Doxygen | 1.8.13              | 1.8.13           |
| MPI     | 2.1, 3.1            | 3.1              |

### Install using CMake: the preferred way

The preferred way to install Tasmanian is to use the included CMake build script, which requires CMake version 3.19,
ROCm capabilities require CMake 3.21.

* The commands for an out-of-source CMake build:
```
  mkdir Build
  cd Build
  cmake <options> <path-to-Tasmanian-source>
  make
  make test
  make install
  make test_install
```

* Standard CMake options are accepted:
```
  -D CMAKE_INSTALL_PREFIX:PATH=<install-prefix> (install folder for the make install command)
  -D CMAKE_BUILD_TYPE:STRING=<Debug/Release>    (set debug flags or default optimization flags)
  -D BUILD_SHARED_LIBS:BOOL=<ON/OFF>            (pick shared or static libs)
  -D CMAKE_CXX_COMPILER:PATH=<path>             (specify the C++ compiler, or HIP/ROCm compiler)
  -D CMAKE_CUDA_COMPILER:PATH=<path>            (specify the CUDA nvcc compiler)
  -D CMAKE_CXX_FLAGS:STING=<flags>              (set additional flags)
```

* List of Tasmanian specific CMake options (all default to **OFF**):
```
  -D Tasmanian_ENABLE_OPENMP:BOOL=<ON/OFF>      (recommended)
  -D Tasmanian_ENABLE_BLAS:BOOL=<ON/OFF>        (recommended)
  -D Tasmanian_ENABLE_PYTHON:BOOL=<ON/OFF>      (recommended)
  -D Tasmanian_ENABLE_RECOMMENDED:BOOL=<ON/OFF> (enable the above and the -O3 flag)
  -D Tasmanian_ENABLE_CUDA:BOOL=<ON/OFF>        (stable)
  -D Tasmanian_ENABLE_HIP:BOOL=<ON/OFF>         (stable)
  -D Tasmanian_ENABLE_DPCPP:BOOL=<ON/OFF>       (mostly stable)
  -D Tasmanian_ENABLE_MAGMA:BOOL=<ON/OFF>       (stable)
  -D Tasmanian_MATLAB_WORK_FOLDER:PATH=""       (stable)
  -D Tasmanian_ENABLE_DOXYGEN:BOOL=<ON/OFF>     (stable)
  -D Tasmanian_ENABLE_FORTRAN:BOOL=<ON/OFF>     (mostly stable)
  -D Tasmanian_ENABLE_MPI:BOOL=<ON/OFF>         (stable)
```

* Acceleration options:
    * OpenMP allows Tasmanian to use more than one CPU core, which greatly increases the performance
    * Basic Linear Algebra Subroutines (BLAS) is a standard with many implementations,
      e.g., [https://www.openblas.net/](https://www.openblas.net/); optimized BLAS improves the
      performance when using evaluate commands on grids with many points or working with models with many outputs
    * Linear Algebra PACKage (LAPACK) is a set of advanced solver, eigen-solvers and decomposition methods
      that build upon BLAS and is usually included in the same package, e.g., OpenBLAS, MKL and ESSL;
      within Tasmanian, the name BLAS in CMake or run-time options indicate the dependence and usage of both BLAS and LAPACK
    * CUDA is a C++ language extension that allows Tasmanian to leverage the computing power of Nvidia GPU devices,
      which greatly enhances the performance of `evaluateFast()` and `evaluateBatch()` and a few other calls
    * HIP/ROCm is very similar to CUDA but uses AMD GPU devices instead, Tasmanian supports a HIP backend
    * SYCL/DPC++/OneAPI is the Intel alternative to CUDA and HIP, Tasmanian supports a SYCL/DPC++ backend
    * Matrix Algebra on GPU and Multicore Architectures (MAGMA) is a library for GPU accelerated linear
      algebra developed at the University of Tennessee at Knoxville
    * MPI allows the use of distributed memory in Bayesian inference, parallel surrogate model construction, and send/receive grid through an MPI communicators
* The **Tasmanian_ENABLE_RECOMMENDED** option searches for OpenMP, BLAS, and Python, enables the options (if possible) and also sets the `-O3` flag
* **Tasmanian_ENABLE_FORTRAN** enables the Fortran 2003 interface
    * Fortran 90 interface is deprecated, but can be enabled with **-DTasmanian_ENABLE_FORTRAN90=ON**
* Additional interfaces are available, beyond the default C/C++ library and the command line tools:
    * Python and Fortran require appropriate interpreter and compiler
    * The MATLAB/Octave interface requires a work-folder with read/write permission for temporary files.
      The interface is enabled by setting **Tasmanian_MATLAB_WORK_FOLDER** to a valid read/write location.
* The Doxygen option will build the HTML documentation

* Options to adjust the testing environment: by default Tasmanian testing will use the system provided OpenMP parameters and run tests on all visible GPU devices; specific number of threads and device can be selected (note that only the testing environment is affected):
```
 -D Tasmanian_TESTS_OMP_NUM_THREADS=<number-of-threads-for-testing> (only used with OpenMP)
 -D Tasmanian_TESTS_GPU_ID=<cuda-device-id-for-testing>             (only used with CUDA)
```

* Additional commands to guide the CMake `find_package()` modules:
```
  -D Python_EXECUTABLE:PATH            (specify the Python interpreter)
  -D CMAKE_CUDA_COMPILER:PATH          (specify the CUDA nvcc compiler)
  -D Tasmanian_ENABLE_ROCM:PATH        (specify the search path for HIP)
  -D CMAKE_Fortran_COMPILER:PATH       (specify the Fortran compiler)
  -D MKLROOT:PATH                      (specify the path to the oneMKL installation)
  -D Tasmanian_MAGMA_ROOT:PATH         (specify the path to the MAGMA installation)
  -D MPI_CXX_COMPILER:PATH=<path>      (specify the MPI compiler wrapper)
  -D MPI_Fortran_COMPILER:PATH=<path>  (needed for MPI with Fortran)
  -D MPIEXEC_EXECUTABLE:PATH=<path>    (needed for MPI testing)
```

* The MAGMA option can be passed without the *Tasmanian_* prefix and if not specified explicitly, it will be read from the OS environment.

* Option for automated download of MAGMA, the package will be downloaded, build, and installed together with Tasmanian in the same prefix.
```
  -D Tasmanian_MAGMA_DOWNLOAD:BOOL=<ON/OFF> (works only on Linux)
```

* The **ROCm** capabilities **NO LONGER** require the CMake CXX compiler to be set to *hipcc*.
    * by default, Tasmanian will search for hip in `/opt/rocm /opt/rocm/hip` following the [ROCm documentation](https://rocmdocs.amd.com/en/latest/Installation_Guide/Using-CMake-with-AMD-ROCm.html)
    * additional search paths can be added using `CMAKE_PREFIX_PATH`
    * the CXX compiler, if not set automatically, should be set to ROCm clang, e.g., `/opt/rocm/llvm/bin/clang++`
    * if `Tasmanian_ENABLE_ROCM` is defined and set to a directory, it will be added to `CMAKE_PREFIX_PATH` to search for `roc::hip`

* The **OneAPI** capabilities require:
    * the compiler has to be set to Intel *ipcx*, Tasmanian will automatically add the *-fsycl* flag (when needed)
    * `Tasmanian_ENABLE_BLAS` is set to **ON**
    * BLAS is set to the CPU version of MKL, e.g., using `BLAS_LIBRARIES` or `BLA_VENDOR` or `MKLROOT` environment variable

* Alternatives allowing to directly specify libraries and bypass `find_package()` altogether:
```
  -D BLAS_LIBRARIES
  -D LAPACK_LIBRARIES
```

* Extra options are available in case CMake fails to find a required dependency, e.g., `find_package()` sometimes fails to acknowledge that the ACML implementation of BLAS depends on both `libgfortran` and `libgomp`; the manual options below should not be necessary in the majority of cases:
```
  -D Tasmanian_EXTRA_LIBRARIES:LIST     (add more libraries as dependencies)
  -D Tasmanian_EXTRA_INCLUDE_DIRS:PATH  (add more include paths to search for headers)
  -D Tasmanian_EXTRA_LINK_DIRS:PATH     (appends more link paths to search for libraries)
```

* Options helpful to Tasmanian developers:
```
  -D DOXYGEN_INTERNAL_DOCS=YES  (include the documentation of the Tasmanian internals)
```

### Install with (basic) GNU Make

The core capabilities of Tasmanian can be build with a few simple GNU Make commands.
The basic build engine is useful for quick testing and exploring Tasmanian, or
if CMake is unavailable or unwanted.
Acceleration options other than OpenMP are not supported in the basic mode,
CMake is so common these days that it should be the default way to install.

* Using GNU Make with `g++` and optionally `/usr/bin/env python`
```
  make
  make test     (will fail if /usr/bin/env python is missing the numpy or ctypes modules)
  make matlab   (optional: sets matlab work folder to ./tsgMatlabWorkFolder/)
  make examples
  make clean
```
In the basic mode, the source folder will become the installation folder, i.e.,
the libraries, executables and Python modules will be build in the source folder
and the headers will be copied to the `include` sub-folder.

### Install with Python Pip

Tasmanian is included in the Python Pip index: [https://pypi.org/project/Tasmanian/](https://pypi.org/project/Tasmanian/)
```
  python3 -m pip install Tasmanian --user   (user installation)
  python3 -m pip install Tasmanian          (virtual env installation)
```

The Tasmanian module is not a regular Python-only project but a wrapper around the C++ libraries, note the following:
* The compiled binaries are cached by pip and `--no-cache-dir` option must be used to change the options.
* Pip versions prior to 1.10 require that dependencies are installed manually.
* Only user installations are supported, installation for all users is possible with CMake but not Pip.
* Python virtual environments are supported, as well as Linux, Mac and Windows operating systems.
* The Pip installer will accept Tasmanian options specified in the environment variables:
```
Environment Option will translate to               CMake Options
export Tasmanian_ENABLE_BLAS=<blas-lapack-libs>    -D Tasmanian_ENABLE_BLAS=ON
                                                   -D BLAS_LIBRARIES=<blas-lapack-libs>
                                                   -D LAPACK_LIBRARIES=<blas-lapack-libs>
export Tasmanian_ENABLE_CUDA=<cuda-nvcc>           -D Tasmanian_ENABLE_CUDA=ON
                                                   -D CMAKE_CUDA_COMPILER=<cuda-nvcc>
export Tasmanian_ENABLE_DPCPP=<dpcpp>              -D Tasmanian_ENABLE_DPCPP=ON
                                                   -D CMAKE_CXX_COMPILER=<dpcpp>
export Tasmanian_ENABLE_MAGMA=<magma-root>         -D Tasmanian_ENABLE_MAGMA=ON
                                                   -D Tasmanian_MAGMA_ROOT_DIR=<magma-root>
export Tasmanian_ENABLE_MPI=<mpicxx>               -D Tasmanian_ENABLE_MPI=ON
                                                   -D MPI_CXX_COMPILER=<mpicxx>
export Tasmanian_MATLAB_WORK_FOLDER=<path>         -D Tasmanian_MATLAB_WORK_FOLDER=<path>
```
* CMake also accepts several standard environment variables:
```
Environment Option          Example values
export BLA_VENDOR           OpenBLAS or MKL
export MKLROOT              /opt/intel/oneapi/mkl/latest
```
* note that the Tasmanian specific environment variables work only with Python Pip
* Example virtual install looks like this:
```
python3 -m venv tasmanian_virtual_env                   # create a virtual environment
source ./tasmanian_virtual_env/bin/activate             # activate the virtual environment
export Tasmanian_ENABLE_CUDA=/usr/local/cuda/bin/nvcc   # specify the CUDA compiler
python -m pip install Tasmanian --no-cache-dir          # will install Tasmanian with CUDA
python -m Tasmanian                                     # print the install log
tasgrid -v                                              # print the available CUDA devices
```

Additional notes:
* under MS Windows the environment variables can be set from Advanced System Settings,
  but the paths should use Linux style of back-slashes, e.g., `C:/Program Files/CUDA/bin/nvcc.exe`
* scikit build sometimes lags in support for the latest MS Visual Studio, regular CMake install works fine
* under OSX some users have reported segfaults when using a pip install, the problem does not
  appear when using CMake, the issue is under investigation


### Install with Spack

Tasmanian is also included in Spack: [https://spack.io/](https://spack.io/)
```
 spack install tasmanian@8.0+openmp+blas+cuda+magma+python+fortran
```

### Install on MS Windows platform

Tasmanian has been tested with MS Visual Studio 2019.

* First use the CMake GUI to set the folders and options
* Then use the command prompt (`cmd.exe`) to enter the build folder
```
  cd <cmake-build-folder>
  cmake --build . --config Release
  ctest -C Release
  cmake --build . --config Release --target install
```
* Both Debug and Release are supported config modes, but do not use them simultaneously,
  pick only one Release or Debug.

### Install folder structure

Tasmanian follows standard Linux conventions, the install path could
be set to `/usr/local/`, although it is recommended to install in
a location inside the user home folder to avoid potential system-wide conflicts.

* Install folder structure:
```
  <install-path>/bin/                       (tagrid executable tool)
  <install-path>/lib/                       (shared or static libraries)
  <install-path>/lib/Tasmanian/             (cmake package-config files)
  <install-path>/lib/pythonX.Y/             (python module)
  <install-path>/include/                   (headers .h and .hpp, and Fortran .mod)
  <install-path>/share/Tasmanian            (bash env scripts, install log, table)
  <install-path>/share/Tasmanian/examples/  (reference examples)
  <install-path>/share/Tasmanian/matlab/    (matlab scripts)
  <install-path>/share/Tasmanian/python/    (copy of <install-path>/lib/pythonX.Y)
  <install-path>/share/Tasmanian/testing/   (a CMake script for post-install testing)
```

Additional notes:
* A summary of all compile options is stored in
```
  <install-path>/share/Tasmanian/Tasmanian.log
```
* The log file can be displayed with wither command:
```
  <install-path>/bin/tasgrid -log  (using the command line tool)
  python -m Tasmanian              (only if Python is enabled)
```
* The executable and library paths, as well as the Python path can be set in `bash` by sourcing
```
  source <install-path>/share/Tasmanian/TasmanianENVsetup.sh
```
* When using pip, the python module location may be different depending on the OS and the Python virtual environment.
* The Python module is version independent, i.e., the files work with all tested Python versions; the `share/Tasmanian/python` copy allows for easy and version independent access to the Tasmanian modules.
* Under MS Windows the shared libraries (e.g., the .dll files) are installed in `bin`.

### Linking to Tasmanian: CMake Package Config

Tasmanian will install CMake package-config files in `<install-path>/lib/Tasmanian`, the files will contain all necessary information to import the Tasmanian targets into other CMake projects using the CMake command:
```
  find_package(Tasmanian 7.0 PATHS "<install-path>")
```
See the included `CMakeLists.txt` in `<install-path>/share/Tasmanian/examples`.

Note that the `PATHS` do not have to be specified explicitly
if the `TasmanianENVsetup.sh` is sourced or if the `Tasmanian_ROOT` environment variable is set.
The correct `find_package()` command is displayed in the log (see the previous section).

The imported targets will be named:
```
  Tasmanian::Tasmanian  (links to the C++ libraries)
  Tasmanian::tasgrid    (imported executable pointing to the command line tool)
  Tasmanian::Fortran    (links to the C++ libraries and the Fortran wrappers)
```

In addition, the following variables will be set:
```
  Tasmanian_PYTHONPATH         (path to the Python module, if Python was enabled)
  Tasmanian_MATLAB_WORK_FOLDER (path to the MATLAB work folder, if set during build)
  Tasmanian_MATLABPATH         (path to the MATLAB scripts, if MATLAB was enabled)
  Tasmanian_<component>_FOUND  (set to ON for each available component)
```
The possible components are:
```
  SHARED STATIC OPENMP BLAS CUDA HIP DPCPP MAGMA MPI PYTHON MATLAB FORTRAN
```
The modules correspond to shared and static libraries and the cmake options used during build.

All available components will be included even if the component is not explicitly requested.
Requesting components can help catch errors early in the build process
and/or print useful log messages. For example:
```
  find_package(Tasmanian 7.0 REQUIRED SHARED PYTHON CUDA OPTIONAL_COMPONENTS OPENMP)
```
In the above example:
* an error will be generated if Tasmanian was build with static libraries, no CUDA or no Python support
* a status message will report whether Tasmanian was build with OpenMP support
* requesting incompatible components will always fail, e.g., CUDA with ROCM or SHARED and STATIC

### Known Issues

Several known issues and work-around fixes:
* the automated MAGMA download reports the wrong GPU architecture enabled
    * ignore the MAGMA message the build correctly respects `CMAKE_CUDA_ARCHITECTURES`
* The addon tests sometime fail due to thread scheduling
    * The overhead associated with thread scheduling is much larger than the simple test models used,
    which leads to unrealistically large fluctuations in sample run-time, which in turn leads to
    randomness in the results, most notably on machines with few cpu cores.
    * Rerun the tests and/or installation to see if the problem is persistent
* Mixing the GCC and Clang compilers and linkers sometimes fails with an error about the architecture
    * use shared libraries only, i.e., `-D BUILD_SHARED_LIBS=ON` in CMake
* Some older versions of the PGI compiler fails when using optimization `-O2`
    * use the `-O1` instead, or the newest version of the compiler
* XL compiler with OpenMP segfaults if the OMP_NUM_THREADS is not set correctly
* XL Fortran needs to use the default halt level which is unnecessarily modified by CMake
    * overwrite the CMake flags using `-DCMAKE_Fortran_FLAGS=-qthreaded`
    * see [the github discussion for details](https://github.com/xsdk-project/xsdk-issues/issues/94#issuecomment-916890894)
* Older versions of CUDA do not work with newer versions of some compilers, e.g., `gcc`
    * consult the CUDA manual for a list of acceptable compilers
