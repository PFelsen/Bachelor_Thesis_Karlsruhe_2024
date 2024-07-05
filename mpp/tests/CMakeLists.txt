add_library(LIB_TEST STATIC TestEnvironment.cpp)
if (BUILD_UQ)
    target_link_libraries(LIB_TEST LIB_UQ ${GTEST_LIB})
else ()
    target_link_libraries(LIB_TEST LIB_APP ${GTEST_LIB})
endif ()

# Test macro to add mpp tests
macro(add_mpp_test testpath lib)
    # To remove pathnames infront of testname
    string(FIND ${testpath} "/" pos REVERSE)
    set(optionalSources ${ARGN})
    MATH(EXPR pos "${pos}+1")
    string(SUBSTRING ${testpath} ${pos} -1 testname)
    add_executable(${testname} ${testpath}.cpp ${optionalSources})
    target_link_libraries(${testname} ${lib} LIB_TEST)
    add_test(${testname} ${testname})
endmacro()

file(WRITE ${PROJECT_BINARY_DIR}/mpi_tests.txt "")
macro(add_mpi_test testpath lib)
    # To remove pathnames infront of testname
    string(FIND ${testpath} "/" pos REVERSE)
    MATH(EXPR pos "${pos}+1")
    string(SUBSTRING ${testpath} ${pos} -1 testname)
    string(APPEND testname "MPI")
    add_executable(${testname} ${testpath}.cpp)
    target_link_libraries(${testname} ${lib} LIB_TEST)
    # Creates file registering mpi tests
    file(APPEND ${PROJECT_BINARY_DIR}/mpi_tests.txt
            "${CMAKE_CURRENT_BINARY_DIR}/${testname}\n")
endmacro()

file(WRITE ${PROJECT_BINARY_DIR}/mpp_bench.txt "")
macro(add_mpp_bench testpath lib)
    # To remove pathnames infront of testname
    string(FIND ${testpath} "/" pos REVERSE)
    MATH(EXPR pos "${pos}+1")
    string(SUBSTRING ${testpath} ${pos} -1 testname)
    string(APPEND testname "BENCH")
    add_executable(${testname} ${testpath}.cpp)
    target_link_libraries(${testname} ${lib} LIB_TEST)
    file(APPEND ${PROJECT_BINARY_DIR}/mpp_bench.txt
            "${CMAKE_CURRENT_BINARY_DIR}/${testname}\n")
endmacro()

# Adds warning to ctest if space dimension is lower than 3, i.e, not all tests are running
MACRO(add_ctest_warning path)
    if (SPACE_DIM MATCHES "^[1-2]$")
        file(WRITE ${path}/CTestCustom.cmake
            "set(CTEST_CUSTOM_PRE_TEST \"echo \nWARNING: Probably not all test run due to lower space dimension...\n\")\n
            set(CTEST_CUSTOM_POST_TEST \"echo \nWARNING: Probably not all test ran due to lower space dimension!\n\")")
    else ()
        file(WRITE ${path}/CTestCustom.cmake "")
    endif ()
ENDMACRO()

macro(add_test_subdirectory testpath)
    add_subdirectory(${testpath})
    add_ctest_warning(${CMAKE_CURRENT_BINARY_DIR}/${testpath})
endmacro()

if (BUILD_TESTS)
    add_ctest_warning(${CMAKE_CURRENT_BINARY_DIR})
    add_ctest_warning(${CMAKE_CURRENT_BINARY_DIR}/..)
    add_test_subdirectory(test0_basic)
    add_test_subdirectory(test1_math)
    add_test_subdirectory(test2_mesh)
    add_test_subdirectory(test3_disc)
    add_test_subdirectory(test4_fem)
    add_test_subdirectory(test5_solve)
endif ()
