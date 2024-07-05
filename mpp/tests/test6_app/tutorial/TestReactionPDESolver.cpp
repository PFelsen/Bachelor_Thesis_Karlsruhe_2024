#include "TestReactionPDESolver.hpp"

/*
 * Todo:
 *  - actually test values!
 *  - test cases for dg and artificial diffusion
 *
 */


TEST_PROBLEMS(TestReactionProblems)

TEST_PROBLEMS(TestReactionProblemsFullCommSplit)

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithoutDefaultConfig()
                     .WithConfPath(std::string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(string(ProjectMppDir) + "/conf/geo/")
                     .WithScreenLogging()
                     .WithFileLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}