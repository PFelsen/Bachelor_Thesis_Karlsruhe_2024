//
// Created by lstengel on 22.07.22.
//

#include "IVectorValuedAssemble.hpp"
#include "TestEnvironment.hpp"


std::vector<std::string> ellipticProblems{"P1Test2D", "P2Test2D"};

// Test Parameters are degree, level, problemName
class VectorValuedEllipticTest : public TestWithParam<std::tuple<int, int, std::string>> {
protected:
  std::unique_ptr<IVectorValuedProblem> problem;
  std::unique_ptr<IVectorValuedAssemble> assemble;

  std::unique_ptr<Vector> u;

  VectorValuedEllipticTest(std::string model) :
      problem(CreateVectorValuedProblemUnique(std::get<2>(GetParam()))) {
    problem->CreateMeshes(MeshesCreator().WithLevel(std::get<1>(GetParam())).WithPLevel(0));
    assemble = CreateVectorValuedAssembleUnique(*problem, PDESolverConfig()
                                                              .WithModel(model)
                                                              .WithDegree(std::get<0>(GetParam()))
                                                              .WithPLevel(0)
                                                              .WithLevel(std::get<1>(GetParam())));

    u = std::make_unique<Vector>(0.0, assemble->GetSharedDisc());
    assemble->Initialize(*u);
    assemble->PrintInfo();
    Config::PrintInfo();
  }
};

/**************************************************
 *  Lagrange Consistency
 **************************************************/


class VectorValuedLagrangeTest : public VectorValuedEllipticTest {
protected:
  VectorValuedLagrangeTest() : VectorValuedEllipticTest("VectorValuedLagrange") {}
};

class VectorValuedDGTest : public VectorValuedEllipticTest {
protected:
  VectorValuedDGTest() : VectorValuedEllipticTest("DGVectorValuedAssemble") {}
};

TEST_P(VectorValuedLagrangeTest, AssembleConsistency) {
  EXPECT_TRUE(isAssembleConsistent(*assemble, *u));
}

INSTANTIATE_TEST_SUITE_P(Degree_1_Consistency, VectorValuedLagrangeTest,
                         Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

INSTANTIATE_TEST_SUITE_P(Degree_2_Consistency, VectorValuedLagrangeTest,
                         Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

TEST_P(VectorValuedDGTest, AssembleConsistency) {
  EXPECT_TRUE(isAssembleConsistent(*assemble, *u));
}

INSTANTIATE_TEST_SUITE_P(Degree_1_Consistency, VectorValuedDGTest,
                         Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

INSTANTIATE_TEST_SUITE_P(Degree_2_Consistency, VectorValuedDGTest,
                         Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

/**************************************************
 *  EG Consistency
 **************************************************/
class VectorValuedEGTest : public VectorValuedEllipticTest {
protected:
  VectorValuedEGTest() : VectorValuedEllipticTest("EGVectorValuedAssemble") {}
};

TEST_P(VectorValuedEGTest, AssembleConsistency) {
  EXPECT_TRUE(isAssembleConsistent(*assemble, *u));
}

INSTANTIATE_TEST_SUITE_P(Degree_1_Consistency, VectorValuedEGTest,
                         Combine(ValuesIn({1}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

INSTANTIATE_TEST_SUITE_P(Degree_2_Consistency, VectorValuedEGTest,
                         Combine(ValuesIn({2}), ValuesIn({0, 2}), ValuesIn(ellipticProblems)));

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfPath(std::string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(std::string(ProjectMppDir) + "/conf/geo/")
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
