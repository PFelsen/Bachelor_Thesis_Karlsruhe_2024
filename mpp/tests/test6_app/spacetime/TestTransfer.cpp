#include <map>
#include <memory>
#include <STDGViscoAcousticAssemble.hpp>
#include "GMRES.hpp"
#include "SpaceTime.hpp"
#include "TestEnvironment.hpp"

class TestTransfer : public Test {
protected:
  void SetUp() override {}

  void TearDown() override {}

  void TestRestrictSolveProlongate(DegreePair degree, LevelPair clevel, LevelPair flevel,
                                   const std::string &problem, const std::string &ttype) {
    mout << DOUT(clevel) << DOUT(flevel) << endl;
    std::unique_ptr<Meshes> meshes =
        MeshesCreator().WithPLevel(0).WithLevel(4).WithMeshName("SpaceTimeSquare").CreateUnique();
    STDGViscoAcousticAssemble assemble(*meshes, degree, problem);


    Vector rhs_F(0.0, assemble.GetSharedDisc(), flevel);
    Vector rhs_C(0.0, assemble.GetSharedDisc(), clevel);
    rhs_F.Clear();
    rhs_C.Clear();

    Matrix M_F(rhs_F);
    Matrix M_C(rhs_C);

    assemble.System(M_F, rhs_F);
    assemble.System(M_C, rhs_C);

    std::unique_ptr<SpaceTimeTransfer> transfer = GetSpaceTimeTransfer(assemble.GetDisc(), ttype);

    auto S_F = GetLinearSolver("GMRES", GetPC("PointBlockGaussSeidel"));
    auto S_C = GetLinearSolver("GMRES", GetPC("PointBlockGaussSeidel"));

    (*S_F)(M_F);
    (*S_C)(M_C);

    Vector sol_F = (*S_F) * rhs_F;
    Vector sol_C = (*S_C) * rhs_C;

    double l2error_F = assemble.L2Error(sol_F);
    double l2error_C = assemble.L2Error(sol_C);


    Vector rhs_C_restricted(0.0, rhs_C);
    transfer->Restrict(rhs_F, rhs_C_restricted);
    Vector sol_C_restricted = (*S_C) * rhs_C_restricted;
    Vector rhs_C_diff = rhs_C - rhs_C_restricted;

    Vector one(1.0, rhs_C_diff);


    Vector sol_F_prolongated(0.0, rhs_F);
    transfer->Prolongate(sol_F_prolongated, sol_C);

    double scalarproductOfDifference = assemble.L2ScalarProduct(one, rhs_C_diff);
    double l2error_C_restricted = assemble.L2Error(sol_C_restricted);
    double l2error_F_prolongated = assemble.L2Error(sol_F_prolongated);

    mout << DOUT(sol_F) << endl;
    mout << DOUT(sol_C) << endl;
    mout << DOUT(sol_F_prolongated) << endl;


    mout << DOUT(l2error_F_prolongated) << endl;
    mout << DOUT(l2error_C_restricted) << endl;


    EXPECT_NEAR(l2error_F, 0.0, 1e-10);
    EXPECT_NEAR(l2error_C, 0.0, 1e-10);

    EXPECT_NEAR(scalarproductOfDifference, 0.0, 1e-12);
    EXPECT_NEAR(l2error_F_prolongated, 0.0, 1e-12);
  }

  void TestAdjoint(DegreePair degree, LevelPair clevel, LevelPair flevel, const std::string &pname,
                   const std::string &tname) {
    mout << DOUT(clevel) << DOUT(flevel) << endl;
    std::unique_ptr<Meshes> meshes =
        MeshesCreator().WithPLevel(0).WithLevel(4).WithMeshName("SpaceTimeSquare").CreateUnique();
    STDGViscoAcousticAssemble assemble(*meshes, degree, pname);
    Vector sol_C(0.0, assemble.GetSharedDisc(), clevel);
    assemble.get_exact_solution(sol_C);
    Vector sol_F(0.0, assemble.GetSharedDisc(), flevel);
    assemble.get_exact_solution(sol_F);
    Vector sol_C_restricted(0.0, sol_C);
    Vector sol_F_prolongated(0.0, sol_F);
    std::unique_ptr<SpaceTimeTransfer> transfer = GetSpaceTimeTransfer(assemble.GetDisc(), tname);
    transfer->Prolongate(sol_F_prolongated, sol_C);
    transfer->Restrict(sol_F, sol_C_restricted);


    mout << DOUT(sol_F) << endl;
    mout << DOUT(sol_C_restricted) << endl;

    mout << DOUT(sol_F_prolongated) << endl;
    mout << DOUT(sol_C_restricted) << endl;


    double sp1 = sol_C * sol_C_restricted;
    double sp2 = sol_F * sol_F_prolongated;
    EXPECT_NEAR(sp1, sp2, 1e-10);


    double sp11 = assemble.L2ScalarProduct(sol_C, sol_C_restricted);
    double sp22 = assemble.L2ScalarProduct(sol_F, sol_F_prolongated);
  }

  void TestProlongatedSolution(DegreePair degree, LevelPair clevel, LevelPair flevel,
                               const std::string &pname, const std::string &tname) {
    mout << DOUT(clevel) << DOUT(flevel) << endl;
    std::unique_ptr<Meshes> meshes =
        MeshesCreator().WithPLevel(0).WithLevel(4).WithMeshName("SpaceTimeSquare").CreateUnique();
    STDGViscoAcousticAssemble assemble(*meshes, degree, pname);
    Vector rhs_C(0.0, assemble.GetSharedDisc(), clevel);
    Matrix M_C(rhs_C);
    assemble.System(M_C, rhs_C);
    Vector sol_C(0.0, rhs_C);
    assemble.get_exact_solution(sol_C);
    Vector res_C = M_C * sol_C;
    res_C -= rhs_C;

    std::unique_ptr<SpaceTimeTransfer> transfer = GetSpaceTimeTransfer(assemble.GetDisc(), tname);

    Vector rhs_F(0.0, assemble.GetSharedDisc(), flevel);
    Matrix M_F(rhs_F);
    assemble.System(M_F, rhs_F);
    Vector sol_F(0.0, rhs_F);
    transfer->Prolongate(sol_F, sol_C);
    Vector res_F = M_F * sol_F;
    res_F -= rhs_F;

    Vector exact_sol_F(0.0, sol_F);
    assemble.get_exact_solution(exact_sol_F);

    Vector diff = exact_sol_F - sol_F;


    double l2_residual_C = assemble.L2Norm(res_C);
    double l2_residual_F = assemble.L2Norm(res_F);

    double l2_error_C = assemble.L2Error(sol_C);
    double l2_error_F = assemble.L2Error(sol_F);

    double l2_norm_C = assemble.L2Norm(sol_C);
    double l2_norm_F = assemble.L2Norm(sol_F);


    EXPECT_NEAR(l2_error_F, 0.0, 1e-12);
    EXPECT_NEAR(l2_error_C, 0.0, 1e-12);

    EXPECT_NEAR(l2_residual_C, 0.0, 1e-12);
    EXPECT_NEAR(l2_residual_F, 0.0, 1e-12);

    mout << "L2Norm of solution                        : " << l2_norm_C << endl;
    mout << "L2Norm of prolonged solution              : " << l2_norm_F << endl;
    mout << "L2Norm Error of solution                  : " << l2_error_C << endl;
    mout << "L2Norm Error of prolongated solution      : " << l2_error_F << endl;
    mout << "L2Norm of residual of solution            : " << l2_residual_C << endl;
    mout << "L2Norm of residual of prolongated solution: " << l2_residual_F << endl;
  }

  void TestRestrictedSolution(DegreePair degree, LevelPair clevel, LevelPair flevel,
                              const std::string &pname, const std::string &tname) {
    mout << DOUT(clevel) << DOUT(flevel) << endl;
    std::unique_ptr<Meshes> meshes =
        MeshesCreator().WithPLevel(1).WithLevel(4).WithMeshName("SpaceTimeSquare").CreateUnique();

    STDGViscoAcousticAssemble assemble(*meshes, degree, pname);

    Vector rhs_F(0.0, assemble.GetSharedDisc(), flevel);
    Matrix M_F(rhs_F);
    assemble.System(M_F, rhs_F);
    Vector sol_F(0.0, rhs_F);
    assemble.get_exact_solution(sol_F);
    Vector res_F = M_F * sol_F;
    res_F -= rhs_F;


    std::unique_ptr<SpaceTimeTransfer> transfer = GetSpaceTimeTransfer(assemble.GetDisc(), pname);


    Vector rhs_C(0.0, assemble.GetSharedDisc(), clevel);
    Matrix M_C(rhs_C);
    assemble.System(M_C, rhs_C);
    Vector sol_C(0.0, rhs_C);
    transfer->Restrict(sol_F, sol_C);
    Vector res_C = M_C * sol_C;
    res_C -= rhs_C;


    mout << DOUT(sol_C) << endl;
    mout << DOUT(sol_F) << endl;

    double l2_residual_C = assemble.L2Norm(res_C);
    double l2_residual_F = assemble.L2Norm(res_F);

    double l2_error_C = assemble.L2Error(sol_C);
    double l2_error_F = assemble.L2Error(sol_F);


    EXPECT_NEAR(l2_error_F, 0.0, 1e-15);
    EXPECT_NEAR(l2_error_C, 0.0, 1e-15);

    EXPECT_NEAR(l2_residual_C, 0.0, 1e-16);
    EXPECT_NEAR(l2_residual_F, 0.0, 1e-16);

    mout << "L2Norm Error of solution                  : " << l2_error_C << endl;
    mout << "L2Norm Error of prolongated solution      : " << l2_error_F << endl;
    mout << "L2Norm of residual of solution            : " << l2_residual_C << endl;
    mout << "L2Norm of residual of prolongated solution: " << l2_residual_F << endl;
  }
};

TEST_F(TestTransfer, TestRestrictSolveProlongateProjectionConstantSpace) {
  TestRestrictSolveProlongate({0, 0}, {0, 0}, {1, 0}, "Polynomial2DDegree0", "Projection");
}

TEST_F(TestTransfer, TestRestrictSolveProlongateProjectionConstantTime) {
  TestRestrictSolveProlongate({0, 0}, {0, 0}, {0, 1}, "Polynomial2DDegree0", "Projection");
}

TEST_F(TestTransfer, TestRestrictSolveProlongateProjectionLinearSpace) {
  TestRestrictSolveProlongate({1, 1}, {0, 0}, {1, 0}, "Polynomial2DDegree0", "Projection");
}

TEST_F(TestTransfer, TestRestrictSolveProlongateProjectionLinearTime) {
  TestRestrictSolveProlongate({1, 1}, {0, 0}, {0, 1}, "Polynomial2DDegree0", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionConstantSpace) {
  TestAdjoint({0, 0}, {0, 0}, {1, 0}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionConstantTime) {
  TestAdjoint({0, 0}, {0, 0}, {0, 1}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionLinearSpace) {
  TestAdjoint({1, 1}, {0, 0}, {1, 0}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionLinearTime) {
  TestAdjoint({1, 1}, {0, 0}, {0, 1}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionQuadraticSpace) {
  TestAdjoint({2, 2}, {0, 0}, {1, 0}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionQuadraticTime) {
  TestAdjoint({2, 2}, {0, 0}, {0, 1}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestAdjointProjectionQuniticSpace) {
  TestAdjoint({5, 5}, {0, 0}, {1, 0}, "Polynomial2DDegree5", "Projection");
}

TEST_F(TestTransfer, TestProlongatedSolutionProjectionConstantSpace) {
  TestProlongatedSolution({0, 0}, {0, 0}, {1, 0}, "Polynomial2DDegree0", "Projection");
}

TEST_F(TestTransfer, TestProlongatedSolutionProjectionConstantTime) {
  TestProlongatedSolution({0, 0}, {0, 0}, {0, 1}, "Polynomial2DDegree0", "Projection");
}

TEST_F(TestTransfer, TestProlongatedSolutionProjectionLinearSpace) {
  TestProlongatedSolution({1, 1}, {0, 0}, {1, 0}, "Polynomial2DDegree1", "Projection");
}

TEST_F(TestTransfer, TestProlongatedSolutionSpaceQuadraticProjection) {
  TestProlongatedSolution({2, 2}, {1, 1}, {2, 1}, "Polynomial2DDegree2", "Projection");
}

TEST_F(TestTransfer, TestProlongatedSolutionTimeQuadraticProjection) {
  TestProlongatedSolution({2, 2}, {1, 1}, {1, 2}, "Polynomial2DDegree2", "Projection");
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithScreenLogging()
                     .WithPPM()
                     .WithoutDefaultConfig()
                     .WithConfigEntry("Overlap", "STCellsWithCorners")
                     .WithConfigEntry("Distribution", "deformed_optimized")
                     .WithConfigEntry("Verbose", "1")
                     .WithConfigEntry("ConfigVerbose", "1")
                     .WithConfigEntry("LinearVerbose", "1")
                     .WithConfigEntry("Steps", "1000")
                     .WithConfigEntry("Epsilon", "1e-14")
                     .WithConfigEntry("Reduction", "1e-13")
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithConfigEntry("useL2Projection", "1"))
      .RUN_ALL_MPP_TESTS();
}
