#include "TestProcSets.hpp"

PROCSETS_INTERVAL_TESTS(TestProcSetsInterval);

PROCSETS_INTERVAL_TESTS(TestProcSetsIntervalWithSplit);

PROCSETS_INTERVAL_TESTS(TestProcSetsIntervalWithDoubleSplit);


PROCSETS_SQUARE_TESTS(TestProcSetsSquare);

PROCSETS_SQUARE_TESTS(TestProcSetsSquareWithSplit);

PROCSETS_SQUARE_TESTS(TestProcSetsSquareWithDoubleSplit);


PROCSETS_HEXAHEDRON_TESTS(TestProcSetsHexahedron);

PROCSETS_HEXAHEDRON_TESTS(TestProcSetsHexahedronWithSplit);

PROCSETS_HEXAHEDRON_TESTS(TestProcSetsHexahedronWithDoubleSplit);

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("DistributionVerbose", 0)
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}