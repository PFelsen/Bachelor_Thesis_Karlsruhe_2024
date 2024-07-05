# Usage of the M++ CMakeLists
#---------------------------------------------------------------------------------------#
# It is absolutely necessary to define a project name in your own project
# using the CMake command project(..) at the beginning of your CMakeLists, since
# otherwise the project path will not be defined correctly. Moreover, the file
# CMakeOptionMacros.txt has to be included directly after the project command, i.e.,
# your CMakeList has to start with the following lines
#   cmake_minimum_required(VERSION ?.?)
#   project(NavierStokes)
#   include(mpp/CMakeOptionMacros.txt)
#
# Afterwards, several options for the build process of M++ (see list below) can be
# specified for your project (if not, the default values from below will be used).
# The options can be set using either the command
#   set_option(<OPT_NAME> <OPT_VALUE>)
# or the sub option command
#   set_sub_option(<OPT_NAME> <SUB_OPT_NAME> <SUB_OPT_VALUE>)
# which sets a sub option <SUB_OPT_NAME> depending on the corresponding
# variable <OPT_NAME>. Note that the option CMAKE_BUILD_TYPE is ignored. Use the
# M++ build types instead
#
# Additionally, the options can be defined by an optional flag when executing CMake,
# e.g., run CMake in the following way
#   cmake .. -D<OPT_NAME>=<OPT_VALUE>
#
# Finally, in your project the CMakeLists file of M++ has to be included using the
# CMake command
#   include(mpp/CMakeLists.txt)
#---------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.15)
# If no project directory is set yet, assume that M++ is the project itself

set(CMAKE_CXX_STANDARD_REQUIRED 20)
set(CMAKE_CXX_STANDARD 20)

if (NOT PROJECT_NAME)
    project(mpp)
    include(CMakeOptionMacros.txt)
    if (USE_SPACETIME)
        set_option(SPACE_DIM 2)
        set_option(BUILD_SPACETIME_TESTS ON)
        set_option(AFFINE_LINEAR_TRAFO ON)
    endif ()
    set(PROJECT_MPP_DIR ${PROJECT_SOURCE_DIR})
    string(TOLOWER "" LCPROJECTNAME)
else ()
    set(PROJECT_MPP_DIR "${PROJECT_SOURCE_DIR}/mpp")
    string(TOLOWER ${PROJECT_NAME} LCPROJECTNAME)
endif ()

add_definitions(-DProjectDir="${LCPROJECTNAME}")

message(STATUS "PROJECT_BINARY_DIR= " ${PROJECT_BINARY_DIR})
message(STATUS "PROJECT_SOURCE_DIR= " ${PROJECT_SOURCE_DIR})
message(STATUS "PROJECT_MPP_DIR=    " ${PROJECT_MPP_DIR})

# In the following the options for building M++ are listed
#---------------------------------------------------------------------------------------#

set_option(BUILD_UQ ON) # [ON, OFF]
if (BUILD_UQ)
    set_option(USE_FFTW ON)
    set_option(USE_SPLITTED_COMMS ON)
    set_option(AFFINE_LINEAR_TRAFO ON)
    set_option(BUILD_TUTORIAL ON)
endif ()

set_option(MPP_BUILD_TYPE MppRelease) # [MppDebugFast, MppDebug]
set_option(DEBUG_LEVEL 1) # TODO: Explain each debug level
set_option(CPP_VERSION 20)
set_option(CMAKE_CXX_STANDARD 20)
set_option(CPP_COMPILER_TYPE c++) # [gnu++]
set_option(SUPPRESS_WARNINGS OFF) # [ON]
set_option(BUILD_TESTS OFF) # [ON] -> M++ unit tests
set_option(CONST_GEOMETRIC_TOLERANCE 1e-10)
set_option(CONST_NEAR_ZERO 1e-15)
set_option(CONST_VERY_LARGE 1e30)
set_option(CONST_INFTY 1e100)
set_option(SPACE_DIM 3) # [1, 2]
set_option(USE_SPACETIME OFF) # [ON]
set_option(AFFINE_LINEAR_TRAFO OFF) # [ON]
set_option(USE_SUPERLU ON) # [OFF]
set_option(USE_CXSC OFF) # [ON, OFF]
set_sub_option(USE_CXSC BUILD_IA_TESTS OFF) # [ON]
set_sub_option(USE_CXSC BUILD_CXSC_TOOLBOX OFF) # [ON]
set_sub_option(BUILD_CXSC_TOOLBOX BUILD_CXSC_EXAMPLES OFF) # [ON]
set_sub_option(BUILD_CXSC_TOOLBOX BUILD_IA_VERIFY_QUAD OFF) # [ON]
set_option(USE_FFTW OFF) # [ON]
set_option(USE_DATAMESH OFF) # [ON]
set_option(REPLACE_GEOMETRIES false) # [true] -> Automatically replace .geo files in legacy format
set_option(BUILD_TUTORIAL ON) # [ON, OFF]
set_sub_option(BUILD_TUTORIAL BUILD_TUTORIAL_TESTS OFF) # [ON, OFF]
set_sub_option(BUILD_UQ BUILD_UQ_TESTS OFF) # [ON, OFF]
set_sub_option(BUILD_UQ USE_TASMANIAN ON) # [ON, OFF]
set_sub_option(BUILD_UQ AGGREGATE_FOR_SOLUTION OFF)

# Compile Vtu Converter (vtk library needed)
set_option(CONVERT_VTU OFF)

# Google test / Google benchmark
#---------------------------------------------------------------------------------------#
add_subdirectory(${PROJECT_MPP_DIR}/submodules/googletest)
include_directories(${PROJECT_MPP_DIR}/submodules/googletest/googletest/include)
include_directories(${PROJECT_MPP_DIR}/submodules/googletest/googlemock/include)
set(BENCHMARK_ENABLE_TESTING OFF CACHE STRING "Disable testing of the benchmark library.")
add_subdirectory(${PROJECT_MPP_DIR}/submodules/benchmark)
include_directories(${PROJECT_MPP_DIR}/submodules/benchmark/include)
set(GTEST_LIB gtest gtest_main gmock gmock_main benchmark::benchmark benchmark::benchmark_main)
#---------------------------------------------------------------------------------------#

# nlohmann/json (https://github.com/nlohmann/json)
#---------------------------------------------------------------------------------------#
include_directories(${PROJECT_MPP_DIR}/submodules/json/single_include/nlohmann)
#---------------------------------------------------------------------------------------#

# Tasmanian
#---------------------------------------------------------------------------------------#
if (BUILD_UQ AND USE_TASMANIAN)
    add_definitions(-DUSE_TASMANIAN)
    set(Tasmanian_ENABLE_RECOMMENDED ON)
    set(Tasmanian_ENABLE_MPI ON)
    add_subdirectory(${PROJECT_MPP_DIR}/submodules/TASMANIAN)
    include_directories(${PROJECT_MPP_DIR}/submodules/TASMANIAN/Config)
    include_directories(${PROJECT_MPP_DIR}/submodules/TASMANIAN/Config/AltBuildSystems)
    include_directories(${PROJECT_MPP_DIR}/submodules/TASMANIAN/SparseGrids)
    include_directories(${PROJECT_MPP_DIR}/submodules/TASMANIAN/SparseGrids/Benchmarks)
    include_directories(${PROJECT_MPP_DIR}/submodules/TASMANIAN/DREAM)
endif ()
#---------------------------------------------------------------------------------------#

# SPRNG5
#---------------------------------------------------------------------------------------#
if (BUILD_UQ)
    if (EXISTS "${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_build.txt" AND
            EXISTS "${PROJECT_BINARY_DIR}/submodules/sprng5/lib/libsprng.a")
        message(STATUS "SPRNG5 found")
    else ()

        message(STATUS "(0/2) Copy sprng5 ${PROJECT_BINARY_DIR}/submodules/sprng5/")
        file(COPY ${PROJECT_MPP_DIR}/submodules/sprng5 DESTINATION "${PROJECT_BINARY_DIR}/submodules")

        message(STATUS "(1/2) Configure sprng5 ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_conf.txt")
        set(ENV{CXX} mpiCC)
        execute_process(COMMAND ./configure --with-mpi=yes --with-fortran=no
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/submodules/sprng5
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_conf.txt OUTPUT_QUIET
                ERROR_FILE ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_conf_err.txt ERROR_QUIET)

        execute_process(COMMAND autoreconf -ivf
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/submodules/sprng5
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_reconf.txt OUTPUT_QUIET
                ERROR_FILE ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_reconf_err.txt ERROR_QUIET)

        message(STATUS "(2/2) Building sprng5 ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_build.txt")
        execute_process(COMMAND make -j --directory=${PROJECT_BINARY_DIR}/submodules/sprng5
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/submodules/sprng5
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_build.txt OUTPUT_QUIET
                ERROR_FILE ${PROJECT_BINARY_DIR}/submodules/sprng5/sprng5_build_err.txt ERROR_QUIET)
    endif ()
    link_directories(${PROJECT_BINARY_DIR}/submodules/sprng5/lib)
    include_directories(${PROJECT_BINARY_DIR}/submodules/sprng5/include)
endif ()

#---------------------------------------------------------------------------------------#

# Blas Lapack
#---------------------------------------------------------------------------------------#
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_definitions(-DLAPACK)
#---------------------------------------------------------------------------------------#

# M++ build type
#---------------------------------------------------------------------------------------#
if (NOT ("${MPP_BUILD_TYPE}" STREQUAL "MppRelease" OR
        "${MPP_BUILD_TYPE}" STREQUAL "MppDebugFast" OR
        "${MPP_BUILD_TYPE}" STREQUAL "MppDebug"))
    message(FATAL_ERROR "Unknown M++ build type ${MPP_BUILD_TYPE}")
endif ()
set(CMAKE_BUILD_TYPE ${MPP_BUILD_TYPE})
#---------------------------------------------------------------------------------------#

# C/C++ compiler options (always -O3 since it is only used for external libraries)
#---------------------------------------------------------------------------------------#
set(CMAKE_C_FLAGS_MPPRELEASE "-O3 -Wno-unused-result -fPIC")
set(CMAKE_C_FLAGS_MPPDEBUGFAST "-O3 -Wno-unused-result -fPIC")
set(CMAKE_C_FLAGS_MPPDEBUG "-O3 -Wno-unused-result -fPIC")
set(CMAKE_CXX_FLAGS "-fPIC -g -std=${CPP_COMPILER_TYPE}${CPP_VERSION}")
set(CMAKE_CXX_FLAGS_MPPRELEASE "-O3 -DDebugLevel=0")
set(CMAKE_CXX_FLAGS_MPPDEBUGFAST "-O3 -DDebugLevel=${DEBUG_LEVEL}")
set(CMAKE_CXX_FLAGS_MPPDEBUG "-O0 -DDebugLevel=${DEBUG_LEVEL}")
#---------------------------------------------------------------------------------------#

# Interval arithmetic options
#---------------------------------------------------------------------------------------#
if (USE_CXSC)
    add_definitions(-DBUILD_IA)

    add_subdirectory(${PROJECT_MPP_DIR}/cxsc)
    include_directories(
            ${PROJECT_MPP_DIR}/cxsc/src
            ${PROJECT_MPP_DIR}/cxsc/src/fi_lib
            ${PROJECT_MPP_DIR}/cxsc/src/rts
            ${PROJECT_MPP_DIR}/cxsc/src/asm)

    include_directories(${PROJECT_MPP_DIR}/src/intervalarithmetic)
endif ()

if (BUILD_CXSC_TOOLBOX)
    message(STATUS "Building cxsc toolbox")
    add_subdirectory(${PROJECT_MPP_DIR}/cxsc/CToolbox/Modules)
    include_directories(${PROJECT_MPP_DIR}/cxsc/CToolbox/Modules)
endif ()

if (BUILD_CXSC_EXAMPLES)
    add_subdirectory(${PROJECT_MPP_DIR}/cxsc/examples)
endif ()

if (BUILD_IA_VERIFY_QUAD)
    add_subdirectory(${PROJECT_MPP_DIR}/tests/test3_disc/quadrature/verifyquadrature)
endif ()
#---------------------------------------------------------------------------------------#

# SuperLU
#---------------------------------------------------------------------------------------#
if (USE_SUPERLU)
    set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
    set(BUILD_SHARED_LIBS ON)
    set(XSDK_ENABLE_Fortran OFF)
    option(enable_internal_blaslib OFF)
    option(enable_matlabmex OFF)
    option(enable_tests OFF)

    include_directories(${PROJECT_MPP_DIR}/submodules/superlu/SRC)
    link_directories(${PROJECT_MPP_DIR}/submodules/superlu/SRC)
    add_subdirectory(${PROJECT_MPP_DIR}/submodules/superlu)

    add_definitions(-DUSE_SUPERLU)
endif ()
#---------------------------------------------------------------------------------------#

# Now use MPI compiler
#---------------------------------------------------------------------------------------#
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
#---------------------------------------------------------------------------------------#

# Warning options
#---------------------------------------------------------------------------------------#
add_compile_options(-Werror=return-type)
if (SUPPRESS_WARNINGS)
    add_compile_options(-Wno-deprecated)
    add_compile_options(-Wno-deprecated-declarations)
endif ()
#---------------------------------------------------------------------------------------#

# Constants
#---------------------------------------------------------------------------------------#
add_definitions(-DGeometricTolerance=${CONST_GEOMETRIC_TOLERANCE})
add_definitions(-DnearToZero=${CONST_NEAR_ZERO})
add_definitions(-DVeryLarge=${CONST_VERY_LARGE})
add_definitions(-Dinfty=${CONST_INFTY})
add_definitions(-DReplaceGeometries=${REPLACE_GEOMETRIES})
add_definitions(-DProjectMppDir="${PROJECT_MPP_DIR}")
add_definitions(-DProjectSourceDir="${PROJECT_SOURCE_DIR}")
add_definitions(-DProjectBuildDir="${PROJECT_BINARY_DIR}")
#---------------------------------------------------------------------------------------#

# Space dimension
#---------------------------------------------------------------------------------------#
if (NOT SPACE_DIM MATCHES "^[1-3]$")
    message(FATAL_ERROR "Space dimension must be a number from 1 to 3")
endif ()
add_definitions(-DSpaceDimension=${SPACE_DIM})
#---------------------------------------------------------------------------------------#

# Spacetime options
#---------------------------------------------------------------------------------------#
if (USE_SPACETIME)
    add_definitions(-DUSE_SPACETIME)
    add_definitions(-DTimeDimension=1)
else ()
    add_definitions(-DTimeDimension=0)
endif ()
#---------------------------------------------------------------------------------------#

# Transformation
#---------------------------------------------------------------------------------------#
if (AFFINE_LINEAR_TRAFO)
    add_definitions(-DAFFINE_TRANSFORMATION)
endif ()

if (BUILD_TUTORIAL)
    add_definitions(-DBUILD_TUTORIAL)
endif ()

if (BUILD_UQ)
    add_definitions(-DBUILD_UQ)
endif ()

if (BUILD_UQ AND AGGREGATE_FOR_SOLUTION)
    add_definitions(-DAGGREGATE_FOR_SOLUTION)
endif ()

# DataMesh options
#---------------------------------------------------------------------------------------#
if (USE_DATAMESH)
    add_definitions(-DUSE_DATAMESH)
endif ()
#---------------------------------------------------------------------------------------#

# Manage directory structure in build directory
#---------------------------------------------------------------------------------------#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/vtk)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/gp)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/py)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/geo)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/vtu)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/log)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/json)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/submodules/sprng5)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/usedconf)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/dual)
#---------------------------------------------------------------------------------------#

# Includes a directory and all subdirectories
#---------------------------------------------------------------------------------------#
MACRO(include_subdirectories lib_dir)
    include_directories(${lib_dir})
    file(GLOB_RECURSE file_paths RELATIVE ${lib_dir} ${lib_dir}/*)
    set(dir_list "")
    FOREACH (file_path ${file_paths})
        get_filename_component(dir_path ${file_path} PATH)
        list(APPEND dir_list ${dir_path})
    ENDFOREACH ()
    list(REMOVE_DUPLICATES dir_list)
    FOREACH (dir ${dir_list})
        include_directories(${lib_dir}/${dir})
    ENDFOREACH ()
ENDMACRO()
#---------------------------------------------------------------------------------------#

# Sets all include directories and adds the subdirectory given by subdir_path
#---------------------------------------------------------------------------------------#
MACRO(add_mpp_subdirectory subdir_path)
    include_subdirectories(${subdir_path})
    add_subdirectory(${subdir_path})
ENDMACRO()
#---------------------------------------------------------------------------------------#

# Find TIRCP (for HoreKa)
#---------------------------------------------------------------------------------------#
find_package(PkgConfig QUIET)
pkg_check_modules(PC_TIRPC libtirpc)

find_path(TIRPC_INCLUDE_DIRS
        NAMES netconfig.h
        PATH_SUFFIXES tirpc
        HINTS ${PC_TIRPC_INCLUDE_DIRS}
)

find_library(TIRPC_LIBRARIES
        NAMES tirpc
        HINTS ${PC_TIRPC_LIBRARY_DIRS}
)
set(TIRPC_VERSION ${PC_TIRPC_VERSION})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(TIRPC
        REQUIRED_VARS TIRPC_LIBRARIES TIRPC_INCLUDE_DIRS
        VERSION_VAR TIRPC_VERSION
)

mark_as_advanced(TIRPC_INCLUDE_DIRS TIRPC_LIBRARIES)

if (TIRPC_FOUND)
    include_directories(${TIRPC_INCLUDE_DIRS})
endif ()
#---------------------------------------------------------------------------------------#

# Find ZLIB (for plot compression)
#---------------------------------------------------------------------------------------#
find_package(ZLIB)
if (ZLIB_FOUND)
    add_definitions(-DVTU11_ENABLE_ZLIB)
endif ()
#---------------------------------------------------------------------------------------#

# Subdirectories
#---------------------------------------------------------------------------------------#
add_mpp_subdirectory(${PROJECT_MPP_DIR}/src)
#---------------------------------------------------------------------------------------#

# Executables
#---------------------------------------------------------------------------------------#
if ("${PROJECT_NAME}" STREQUAL "mpp")
    add_executable(M++ ${PROJECT_MPP_DIR}/src/Main.cpp)
    target_link_libraries(M++ MPP_LIBRARIES LIB_APP)

    if (USE_SPACETIME)
        target_link_libraries(M++ LIB_APP)
    endif ()

    if (BUILD_TUTORIAL)
        target_link_libraries(M++ LIB_APP)
    endif ()

    if (BUILD_UQ)
        target_link_libraries(M++ LIB_UQ)
    endif ()

    if (USE_CXSC)
        add_executable(M++Eigenvalue ${PROJECT_MPP_DIR}/src/lib6_app/spectrum/MainEigenvalue.cpp)
        target_link_libraries(M++Eigenvalue MPP_LIBRARIES)
    endif ()
endif ()
#---------------------------------------------------------------------------------------#

# Tests
#---------------------------------------------------------------------------------------#
include(CTest)
enable_testing()

add_mpp_subdirectory(${PROJECT_MPP_DIR}/tests)
if (BUILD_TUTORIAL AND BUILD_TUTORIAL_TESTS)
    add_mpp_subdirectory(${PROJECT_MPP_DIR}/tests/test6_app/tutorial)
endif ()
if (USE_SPACETIME AND BUILD_SPACETIME_TESTS)
    add_subdirectory(${PROJECT_MPP_DIR}/tests/test6_app/spacetime)
endif ()

##Todo: consider refacotring/renaming -> Also idk if anything happens from this alone?

if (BUILD_UQ AND BUILD_UQ_TESTS)
    add_subdirectory(${PROJECT_MPP_DIR}/tests/test7_uq)# add_mpp or add?
endif ()
#---------------------------------------------------------------------------------------#

# Tools
add_mpp_subdirectory(${PROJECT_MPP_DIR}/tools/geometries)


# Print branch names in commit messages
#---------------------------------------------------------------------------------------#
execute_process(COMMAND bash -c ${PROJECT_MPP_DIR}/doc/addCommitMsgHook.sh)
#---------------------------------------------------------------------------------------#

# Copy python and pipeline directory to binary directory
#---------------------------------------------------------------------------------------#
file(COPY ${PROJECT_MPP_DIR}/python DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_MPP_DIR}/pipeline DESTINATION ${PROJECT_BINARY_DIR})
#---------------------------------------------------------------------------------------#

# TODO: Remove this and refactor mppyrun
file(COPY ${PROJECT_MPP_DIR}/python/mppyrun.py DESTINATION ${PROJECT_BINARY_DIR})
