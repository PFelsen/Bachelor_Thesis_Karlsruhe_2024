add_subdirectory(lib0_basic)
add_subdirectory(lib1_math)
add_subdirectory(lib2_mesh)
add_subdirectory(lib3_disc)
add_subdirectory(lib4_fem)
add_subdirectory(lib5_solve)
add_subdirectory(lib6_app)
if (BUILD_UQ)
    add_subdirectory(lib7_uq)
endif()

# TODO REMOVE m++ cpp here and add it to lib0, DEPENDENCY ON LIB4 HAS TO BE REMOVED IN m++.cpp!!!!
add_library(MPP_LIBRARIES SHARED lib0_basic/m++.cpp)
# TODO RENAME THIS: LIB6/7 is final lib now
target_link_libraries(MPP_LIBRARIES LIB_SOLVE)
