#include "TestPollutionProblem.hpp"


INSTANTIATE_TEST_SUITE_P(TestPollutionProblem, TestSpaceTimePollutionProblem,
                         Values(std::pair{ConfigMap{{"Problem", "PollutionSquare500"},
                                                    {"level", "1"},
                                                    {"plevel", "0"},
                                                    {"degree", "1"},
                                                    {"degree_time", "0"},
                                                    {"T", "1.0"},
                                                    {"dt", "0.05"},
                                                    {"ClearDistribution", "0"},
                                                    {"Distribution", "RCB"},
                                                    {"ConfigVerbose", "100"},
                                                    {"LinearVerbose", "1"},
                                                    {"PDESolverVerbose", "1"}
                                                    },
                                          ValueMap{{"L2_Norm", 2.25618}}}));

INSTANTIATE_TEST_SUITE_P(TestPollutionProblem, TestTimeSteppingPollutionProblem,
                         Values(std::pair{ConfigMap{{"Problem", "PollutionSquare500"},
                                                    {"level", "1"},
                                                    {"plevel", "0"},
                                                    {"degree", "1"},
                                                    {"T", "1.0"},
                                                    {"dt", "0.05"},
                                                    {"ClearDistribution", "0"},
                                                    {"Distribution", "RCB"},
                                                    {"rkorder", "-1"}},
                                          ValueMap{{}}}));

TEST_P(TestSpaceTimePollutionProblem, TestRun) { TestRun(); }

TEST_P(TestTimeSteppingPollutionProblem, TestRun) { TestRun(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfPath(std::string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(string(ProjectMppDir) + "/conf/geo/")
                     .WithoutDefaultConfig()
                     .WithScreenLogging()
                     .WithFileLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
