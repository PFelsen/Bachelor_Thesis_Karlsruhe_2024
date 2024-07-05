#include "TestTransportPDESolver.hpp"

#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <PDESolver.hpp>

#include "TestEnvironment.hpp"

/*
 * Todo Open Test Cases:
 *  - All TransportProblems.cpp should either have a test case or be deleted if unused
 *  - Convergence tests
 *  - combination of solvers and problems
 *
 */

namespace {

const std::vector<std::pair<ConfigMap, ValueMap>> configs =
    {std::pair{ConfigMap{{"level", "8"},
                         {"plevel", "8"},
                         {"degree", "2"},
                         {"Problem", "Riemann1D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 1.0}, {"Energy", 3.9999258206882993}}},

     std::pair{ConfigMap{{"level", "8"},
                         {"plevel", "8"},
                         {"degree", "2"},
                         {"Problem", "CosHat1D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 0.19634958}, {"Energy", 0.073631077}}},

     std::pair{ConfigMap{{"level", "8"},
                         {"plevel", "8"},
                         {"degree", "2"},
                         {"Problem", "GaussHat1D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 1.0}, {"Energy", 2.8209477}}},

     std::pair{ConfigMap{{"level", "4"},
                         {"plevel", "4"},
                         {"degree", "0"},
                         {"Problem", "Riemann2D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 0.9999985257627233}, {"Energy", 1.5627957200831057}}},

     std::pair{ConfigMap{{"level", "4"},
                         {"plevel", "4"},
                         {"degree", "2"},
                         {"Problem", "Riemann2D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 1.0}, {"Energy", 3.7100471618094857}}},

     std::pair{ConfigMap{{"level", "4"},
                         {"plevel", "4"},
                         {"degree", "2"},
                         {"Problem", "CosHat2D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 0.036012846}, {"Energy", 0.010416204195390976}}},

     std::pair{ConfigMap{{"level", "4"},
                         {"plevel", "4"},
                         {"degree", "2"},
                         {"Problem", "GaussHat2D"},
                         {"t0", "0.0"},
                         {"T", "0.125"},
                         {"dt", "0.03125"}},
               ValueMap{{"Mass", 0.99958341}, {"Energy", 7.9407467095329629}}},

     std::pair{ConfigMap{{"Problem", "CircleWave2D"},
                         {"level", "3"},
                         {"plevel", "1"},
                         {"degree", "0"},
                         {"dt", "0.005"},
                         {"t0", "0.0"},
                         {"T", "1.0"},
                         {"rkorder", "4"}},
               ValueMap{{"Mass", 5.9038778099093037}, {"Energy", 0.10605307530872354}}},

     std::pair{ConfigMap{{"Problem", "PollutionSquare500"},
                         {"level", "1"},
                         {"plevel", "0"},
                         {"degree", "1"},
                         {"T", "1.6"},
                         {"dt", "0.01"},
                         {"rkorder", "-2"}},
               ValueMap{{"Mass", 0.0070799576204027105}, {"Energy", 0.0021081823020453437}}}};

} // namespace

INSTANTIATE_TEST_SUITE_P(TestTransportPDESolver, TestProblems, ValuesIn(configs));

TEST_P(TestProblems, TestRun) { TestRun(); }

// Skip PollutionSquare500
INSTANTIATE_TEST_SUITE_P(TestTransportPDESolver, TestProblemsCommSplit,
                         ValuesIn(std::begin(configs), std::prev(std::end(configs), 1)));

TEST_P(TestProblemsCommSplit, TestRun) { TestRun(); }

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