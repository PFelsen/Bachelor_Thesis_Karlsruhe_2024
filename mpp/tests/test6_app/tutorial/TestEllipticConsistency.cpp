#include <utility>

#include "IEllipticAssemble.hpp"
#include "PDESolver.hpp"
#include "TestEnvironment.hpp"

std::vector<std::string> ellipticProblems{"Kellogg",        "Rock",      "Laplace1D",
                                          "Laplace2D",      "Divergent", "Discontinuous1D",
                                          "Discontinuous2D"};

// Test Parameters are degree, level, problemName
class EllipticTest : public TestWithParam<std::tuple<int, int, std::string>> {
protected:
  std::unique_ptr<IEllipticProblem> problem;
  std::unique_ptr<IEllipticAssemble> assemble;

  std::unique_ptr<Vector> u;

  explicit EllipticTest(std::string model) {
    problem = CreateEllipticProblemUnique(std::get<2>(GetParam()));
    assemble =
        CreateEllipticAssembleUnique(*problem, PDESolverConfig().WithModel(std::move(model)));
    u = std::make_unique<Vector>(0.0, assemble->GetSharedDisc());
    assemble->Initialize(*u);
  }
};

/**************************************************
 *  Lagrange Consistency
 **************************************************/
class LagrangeEllipticTest : public EllipticTest {
protected:
  LagrangeEllipticTest() : EllipticTest("Lagrange") {}
};

TEST_P(LagrangeEllipticTest, AssembleConsistency) {
  EXPECT_TRUE(isAssembleConsistent(*assemble, *u));
}

INSTANTIATE_TEST_CASE_P

    (Degree_1_Consistency, LagrangeEllipticTest,
     Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));
INSTANTIATE_TEST_CASE_P

    (Degree_2_Consistency, LagrangeEllipticTest,
     Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

/**************************************************
 *  Mixed Consistency
 **************************************************/
class MixedEllipticTest : public EllipticTest {
protected:
  MixedEllipticTest() : EllipticTest("Mixed") {}
};

TEST_P(MixedEllipticTest, AssembleConsistency) { EXPECT_TRUE(isAssembleConsistent(*assemble, *u)); }

INSTANTIATE_TEST_CASE_P

    (Degree_1_Consistency, MixedEllipticTest,
     Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));
INSTANTIATE_TEST_CASE_P

    (Degree_2_Consistency, MixedEllipticTest,
     Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

/**************************************************
 *  Hybrid Consistency
 **************************************************/
class HybridEllipticTest : public EllipticTest {
protected:
  HybridEllipticTest() : EllipticTest("Hybrid") {}
};

TEST_P(HybridEllipticTest, AssembleConsistency) {
  EXPECT_TRUE(isAssembleConsistent(*assemble, *u));
}

INSTANTIATE_TEST_CASE_P

    (Degree_1_Consistency, HybridEllipticTest,
     Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));
INSTANTIATE_TEST_CASE_P

    (Degree_2_Consistency, HybridEllipticTest,
     Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

/**************************************************
 *  DG Consistency
 **************************************************/
class DGEllipticTest : public EllipticTest {
protected:
  DGEllipticTest() : EllipticTest("DG") {}
};

TEST_P(DGEllipticTest, AssembleConsistency) { EXPECT_TRUE(isAssembleConsistent(*assemble, *u)); }

INSTANTIATE_TEST_CASE_P

    (Degree_1_Consistency, DGEllipticTest,
     Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));
INSTANTIATE_TEST_CASE_P

    (Degree_2_Consistency, DGEllipticTest,
     Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

/**************************************************
 *  EG Consistency
 **************************************************/
// TODO: Fix EG Consistency

class EGEllipticTest : public EllipticTest {
protected:
  EGEllipticTest() : EllipticTest("EG") {}
};

TEST_P(EGEllipticTest, AssembleConsistency) { EXPECT_TRUE(isAssembleConsistent(*assemble, *u)); }

INSTANTIATE_TEST_CASE_P

    (Degree_1_Consistency, EGEllipticTest,
     Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));


INSTANTIATE_TEST_CASE_P

    (Degree_2_Consistency, EGEllipticTest,
     Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfPath(std::string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(std::string(ProjectMppDir) + "/conf/geo/")
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
