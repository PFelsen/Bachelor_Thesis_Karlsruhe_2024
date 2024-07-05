#include "TestAcousticPDESolver.hpp"

TEST_PROBLEMS(TestProblems)

TEST_PROBLEMS(TestProblemsFullCommSplit)

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