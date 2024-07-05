#include <map>
#include <memory>
#include "DebuggingTools.hpp"
#include "GMRES.hpp"
#include "SpaceTime.hpp"
#include "TestEnvironment.hpp"

/*
 * Test proposal
 *
 * class TestPolyProblem : public TestWithParam<std::string> {
 *  Sets up Meshes, assemble and solver
 * }
 *
 * class TestConstProblem : public TestPolyProblem
 * class TestLinearProblem : public TestPolyProblem
 * class TestQuadraticProblem : public TestPolyProblem
 *    ...
 * up to max deg of p-adapt
 *
 * TestParam -> ModelString
 *
 * Test() { solve system compare with exact sol or simply check that err = 0.0 }
 *
 * Execute Same Test for all classes with all ModelStrings
 *
 */


class TestQuadratic : public TestWithParam<std::tuple<std::string, std::string>> {
protected:
  std::unique_ptr<Meshes> meshes;

  std::unique_ptr<STAssemble> assemble;

  std::string disc;
  std::string problem;

  // Todo use new main
  void SetUp() override {
    meshes =
        MeshesCreator().WithMeshName("SpaceTimeSquare").WithPLevel(1).WithLevel(1).CreateUnique();

    disc = std::get<0>(GetParam());
    problem = std::get<1>(GetParam());

    assemble = CreateSTAssemble(disc, *meshes, {2, 2}, problem);
  }

  void TearDown() override {}

  void testSolvingAndExactSolution() {
    Vector RHS(0.0, assemble->GetSharedDisc());
    RHS.Clear();
    Matrix M(RHS);
    assemble->System(M, RHS);
    GMRES S(GetPC("PointBlockGaussSeidel"));
    S(M);
    Vector solution = S * RHS;
    double l2error = assemble->L2Error(solution);
    double l2norm = assemble->L2Norm(solution);
    Vector exact_solution(0.0, assemble->GetSharedDisc());
    assemble->get_exact_solution(exact_solution);
    double l2error_exact = assemble->L2Error(exact_solution);
    double l2norm_exact = assemble->L2Norm(exact_solution);
    Vector difference = solution - exact_solution;

    Vector residual = M * solution;
    residual -= RHS;
    double solution_residual = residual.norm();
    residual = M * exact_solution;
    residual -= RHS;
    double exact_solution_residual = residual.norm();

    mout.PrintInfo(disc + " " + problem, 1,
                   std::vector<PrintInfoEntry<double>>{{"L2Norm(solution)", l2norm},
                                                       {"L2Norm(exact_solution)", l2norm_exact},
                                                       {"norm(solution)", norm(solution)},
                                                       {"norm(exact_solution)",
                                                        norm(exact_solution)},
                                                       {"L2Error(solution)", l2error},
                                                       {"L2Error(exact_solution)", l2error_exact},
                                                       {"norm(difference)", norm(difference)},
                                                       {"norm(solution_residual)",
                                                        solution_residual},
                                                       {"norm(exact_solution_residual)",
                                                        exact_solution_residual}});
    ASSERT_LT(l2error, 1e-10);
    ASSERT_LT(l2error_exact, 1e-10);
    ASSERT_LT(solution_residual, 1e-10);
    ASSERT_LT(exact_solution_residual, 1e-10);
  }
};

INSTANTIATE_TEST_SUITE_P(TestQuadratic, TestQuadratic,
                         testing::Combine(Values("STPGViscoAcousticAssemble",
                                                 "STDGViscoAcousticAssemble",
                                                 "STGLGLViscoAcousticAssemble"),
                                          Values("QuadraticST", "QuadraticDirichlet")));

TEST_P(TestQuadratic, SingleCell) { this->testSolvingAndExactSolution(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithoutDefaultConfig()
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Verbose", -1)
                     .WithConfigEntry("ConfigVerbose", 0)
                     .WithConfigEntry("LinearVerbose", -1)
                     .WithConfigEntry("Steps", 1000)
                     .WithConfigEntry("Epsilon", 1e-14)
                     .WithConfigEntry("Reduction", 1e-13)
                     .WithConfigEntry("LinearEpsilon", 1e-14)
                     .WithConfigEntry("LinearReduction", 1e-13)
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
