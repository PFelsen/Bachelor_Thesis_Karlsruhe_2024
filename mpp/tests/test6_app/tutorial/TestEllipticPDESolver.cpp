#include "TestEllipticPDESolver.hpp"


/*
 * Todo Open Test Cases:
 *  - Mixed / Hybrid tests for Polynomial Problems of degree higher than 0,1
 *      -> RTLagrangeElements of higher order are needed
 *  - Develop test problems with triangular cells
 *  - Split up test problems in exact and inexact and merge disc tests
 */

#if (SpaceDimension < 3)
TEST_POLYNOMIALS_2D(TestPolynomialsLagrange)
#endif

#if (SpaceDimension == 3)
TEST_POLYNOMIALS_3D(TestPolynomialsLagrange)
#endif

TEST_PROBLEMS_LAGRANGE(TestProblemsLagrange)

TEST_PROBLEMS_LAGRANGE(TestProblemsLagrangeCommSplit)

#if (SpaceDimension < 3)
TEST_POLYNOMIALS_2D(TestPolynomialsDG)
#endif

#if (SpaceDimension == 3)
TEST_POLYNOMIALS_3D(TestPolynomialsDG)
#endif

TEST_PROBLEMS_DG(TestProblemsDG)

TEST_PROBLEMS_DG(TestProblemsDGCommSplit)

#if (SpaceDimension < 3)
TEST_POLYNOMIALS_2D(TestPolynomialsEG)
#endif

#if (SpaceDimension == 3)
TEST_POLYNOMIALS_3D(TestPolynomialsEG)
#endif

TEST_PROBLEMS_EG(TestProblemsEG)

TEST_PROBLEMS_EG(TestProblemsEGCommSplit)

TEST_PROBLEMS_MIXED(TestProblemsMixed)

TEST_PROBLEMS_MIXED(TestProblemsMixedCommSplit)

TEST_PROBLEMS_MIXED(TestProblemsHybrid)

TEST_PROBLEMS_MIXED(TestProblemsHybridCommSplit)

TEST_SOLVERS(TestLinearSolver)

TEST_SOLVERS(TestGMRES)

TEST_SOLVERS(TestCG)

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithoutDefaultConfig()
                     .WithConfPath(string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(string(ProjectMppDir) + "/conf/geo/")
                     .WithScreenLogging()
                     .WithFileLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}