add_mpp_test(geometry/TestCoarseGeometry LIB_MESH)

add_mpp_test(cells/TestCellConstruction LIB_MESH)
add_mpp_test(cells/TestCellInterface LIB_MESH)


add_mpp_test(mesh/TestCheckOrientation LIB_MESH)
add_mpp_test(mesh/TestMeshGeneration LIB_MESH)
add_mpp_test(mesh/TestMeshCollections LIB_MESH)
add_mpp_test(meshparts/TestMeshPart LIB_MESH)

add_mpi_test(distribution/TestDistribution LIB_MESH)
add_mpi_test(mesh/TestMeshesCreator LIB_MESH)
add_mpi_test(meshparts/TestProcSets LIB_MESH)
add_mpi_test(mesh/TestMeshes LIB_MESH)
add_mpi_test(mesh/TestRefineMesh LIB_MESH)

if (USE_SPACETIME)
    #add_mpi_test(mesh/TestSpaceTimeMesh LIB_MESH)
endif ()

add_mpp_bench(BenchCellCollections LIB_MESH)