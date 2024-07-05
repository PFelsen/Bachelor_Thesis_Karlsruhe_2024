#ifndef SPACETIME_TESTMEMORYLEAK_HPP
#define SPACETIME_TESTMEMORYLEAK_HPP

#include "ErrorEstimator.hpp"
#include "MeshesCreator.hpp"
#include "STDGViscoAcousticAssemble.hpp"
#include "TestEnvironment.hpp"

class TestConformingInterpolation : public Test {
public:
  void SetUp() override {}

  void TearDown() override {}
};

TEST_F(TestConformingInterpolation, TestConstant) {
  auto meshes = MeshesCreator().WithMeshName("ST_QD").WithPLevel(1).WithLevel(2).CreateUnique();

  auto assemble =
      std::make_unique<STDGViscoAcousticAssemble>(*meshes, DegreePair{1, 1}, "Polynomial2DDegree1");

  Vector vec(1.4142, assemble->GetSharedDisc());
  Vector vec2(0.0, assemble->GetSharedDisc());
  assemble->ConformingLagrangeInterpolation(vec, vec2);
  Vector difference = vec - vec2;
  double norm_difference = norm(difference);
  EXPECT_LE(norm_difference, 1e-10);
}

void TestConformingInterpolationForDegree(short degree) {
  auto meshes = MeshesCreator().WithMeshName("ST_QD").WithPLevel(1).WithLevel(2).CreateUnique();
  auto assemble =
      std::make_unique<STDGViscoAcousticAssemble>(*meshes, DegreePair{degree, degree},
                                                  "Polynomial2DDegree" + to_string(degree));
  Vector vec(0.0, assemble->GetSharedDisc());
  assemble->get_exact_solution(vec);
  Vector vec2(0.0, vec);
  assemble->ConformingLagrangeInterpolation(vec, vec2, degree);
  Vector difference = vec - vec2;
  double norm_difference = norm(difference);
  mout << "Error for Degree " << degree << " : " << norm_difference << endl;
  EXPECT_LE(norm_difference, 1e-10);
}

TEST_F(TestConformingInterpolation, TestDegree1) { TestConformingInterpolationForDegree(1); }

TEST_F(TestConformingInterpolation, TestDegree2) { TestConformingInterpolationForDegree(2); }

TEST_F(TestConformingInterpolation, TestDegree3) { TestConformingInterpolationForDegree(3); }

TEST_F(TestConformingInterpolation, TestDegree4) { TestConformingInterpolationForDegree(4); }

TEST_F(TestConformingInterpolation, TestDegree5) { TestConformingInterpolationForDegree(5); }

TEST_F(TestConformingInterpolation, TestErrorOfConformingIterpolationOfRiemannSolution) {
  globalPlotBackendCreator = std::make_unique<DefaultSpaceTimeBackend>;
  auto meshes = MeshesCreator().WithMeshName("ST_QD").WithPLevel(2).WithLevel(3).CreateUnique();

  auto assemble =
      std::make_unique<STDGViscoAcousticAssemble>(*meshes, DegreePair{1, 1}, "DoubleRiemann");

  assemble->GetDisc().GetMeshes().PrintInfo();
  Vector sol_projection(0.0, assemble->GetSharedDisc());
  assemble->get_projected_exact_solution(sol_projection);
  std::unique_ptr<LinearSolver> solver =
      GetLinearSolverUnique("GMRES",
                            std::unique_ptr<Preconditioner>(GetPC("PointBlockGaussSeidel")), "");

  Vector RHS(0.0, sol_projection);
  RHS.Clear();
  Matrix M(sol_projection);
  M = 0;
  assemble->System(M, RHS);
  (*solver)(M);
  L2ReliableResidualErrorEstimator ee(*assemble, *solver);


  Vector sol(0.0, sol_projection);
  sol = (*solver) * RHS;


  Vector conf_sol(0.0, sol_projection);
  assemble->ConformingLagrangeInterpolation(sol, conf_sol);

  Vector conf_proj(0.0, sol_projection);
  assemble->ConformingLagrangeInterpolation(sol_projection, conf_proj);

  Vector eta_sol = ee.EstimateError(M, sol, 0);
  double eta_sol_value = eta_sol.norm();


  Vector conf_err = sol - conf_sol;
  Vector eta = ee.EstimateError(M, conf_proj, 0);
  double eta_value = eta.norm();
  mout << DOUT(eta_value) << endl;
}

#endif // SPACETIME_TESTMEMORYLEAK_HPP
