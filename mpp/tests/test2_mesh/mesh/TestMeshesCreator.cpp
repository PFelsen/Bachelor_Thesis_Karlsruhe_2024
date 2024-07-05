#include "TestMeshesCreator.hpp"

MESHESCREATOR_INTERVAL_TESTS(TestMeshesCreatorInterval);

MESHESCREATOR_INTERVAL_TESTS(TestMeshesCreatorIntervalWithSplit);

MESHESCREATOR_INTERVAL_TESTS(TestMeshesCreatorIntervalWithDoubleSplit);

#if (SpaceDimension >= 2)
MESHESCREATOR_TRIANGLE_TESTS(TestMeshesCreatorTriangle);

MESHESCREATOR_TRIANGLE_TESTS(TestMeshesCreatorTriangleWithSplit);

MESHESCREATOR_TRIANGLE_TESTS(TestMeshesCreatorTriangleWithDoubleSplit);


MESHESCREATOR_SQUARE_TESTS(TestMeshesCreatorSquare);

MESHESCREATOR_SQUARE_TESTS(TestMeshesCreatorSquareWithSplit);

MESHESCREATOR_SQUARE_TESTS(TestMeshesCreatorSquareWithDoubleSplit);


MESHESCREATOR_SQUARETRIANGLE_TESTS(TestMeshesCreatorSquareTriangle);

MESHESCREATOR_SQUARETRIANGLE_TESTS(TestMeshesCreatorSquareTriangleWithSplit);

MESHESCREATOR_SQUARETRIANGLE_TESTS(TestMeshesCreatorSquareTriangleWithDoubleSplit);
#endif

#if (SpaceDimension >= 3)
MESHESCREATOR_TETRAHEDRON_TESTS(TestMeshesCreatorTetrahedron);

MESHESCREATOR_TETRAHEDRON_TESTS(TestMeshesCreatorTetrahedronWithSplit);

MESHESCREATOR_TETRAHEDRON_TESTS(TestMeshesCreatorTetrahedronWithDoubleSplit);


MESHESCREATOR_HEXAHEDRON_TESTS(TestMeshesCreatorHexahedron);

MESHESCREATOR_HEXAHEDRON_TESTS(TestMeshesCreatorHexahedronWithSplit);

MESHESCREATOR_HEXAHEDRON_TESTS(TestMeshesCreatorHexahedronWithDoubleSplit);
#endif


#if (defined USE_SPACETIME && SpaceDimension >= 2)
MESHESCREATOR_SPACETIMESQUARE_TESTS(TestMeshesCreatorSpaceTimeSquare);

MESHESCREATOR_SPACETIMESQUARE_TESTS(TestMeshesCreatorSpaceTimeSquareWithSplit);

MESHESCREATOR_SPACETIMESQUARE_TESTS(TestMeshesCreatorSpaceTimeSquareWithDoubleSplit);
#endif

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("TestMeshesCreatorVerbose", 0)
                     .WithConfigEntry("DistributionVerbose", 0)
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}