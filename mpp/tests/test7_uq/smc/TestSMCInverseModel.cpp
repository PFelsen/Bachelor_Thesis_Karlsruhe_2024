#include "IEllipticAssemble.hpp"
#include "LagrangeEllipticAssemble.hpp"
#include "Newton.hpp"
#include "PDEModel.hpp"
#include "TestEnvironment.hpp"
#include "EllipticProblems.cpp"


TEST(PDEModel, ObservationOperator) {
  std::shared_ptr<RainLaplace2D> problem =
      std::make_shared<RainLaplace2D>(ProposalConfig().WithRandomField("Exp2DKL"));
  problem->GetInput().GetRandomField().GetParameterVector() = RVector(0.0, 1);
  problem->CreateMeshes(MeshesCreator().WithPLevel(1).WithLevel(2));
  std::shared_ptr<IEllipticAssemble> assemble =
      std::make_shared<LagrangeEllipticAssemble>(*problem, 1);
  std::shared_ptr<Newton> newton =
      std::make_shared<Newton>(std::unique_ptr<LinearSolver>(GetLinearSolver()));
  Vector u(0.0, assemble->GetSharedDisc());
  u.SetAccumulateFlag(true);
  newton->operator()(*assemble, u);

  ObservationOperator observation({{"Point", {0.25, 0.25}},
                                   {"Point", {0.50, 0.25}},
                                   {"Point", {0.75, 0.25}},
                                   {"Point", {0.25, 0.50}},
                                   {"Point", {0.50, 0.50}},
                                   {"Point", {0.75, 0.50}},
                                   {"Point", {0.25, 0.75}},
                                   {"Point", {0.50, 0.75}},
                                   {"Point", {0.75, 0.75}}});

  RVector y = observation.GetObservation(u);

  EXPECT_FLOAT_EQ(y[0], -0.25);
  EXPECT_FLOAT_EQ(y[3], -0.50);
  EXPECT_FLOAT_EQ(y[6], -0.75);
}

TEST(PDEModel, EllipticInverseProblem1DSecond) {
  EllipticPDEModel model(
      PDESolverConfig()
          .WithObservations({{"Point", {0.50, 0.50}}})
          .WithModel("Lagrange")
          .WithProblem("RainLaplace2D")
          .WithProposalConfig(ProposalConfig().WithRandomField("Exp2DKL").WithPriorConfig(
              PriorSamplerConfig().WithGenerator("Normal").WithMean(0.0).WithVariance(1.0))));
  model.GetProblem().GetInput().GetRandomField().GetParameterVector() = RVector(0.0, 1);
  RVector u;
  model.Run(u);

  EXPECT_NEAR(u[0], -0.5, 1e-6);
}

TEST(PDEModel, TestModel) {
  TestModel pde_model(
      PDESolverConfig().WithProposalConfig(ProposalConfig().WithRandomField("Simple2DKL")));
  pde_model.GetProblem().GetInput().GetRandomField().GetParameterVector() = RVector(1.0, 1);
  RVector solution(1);
  pde_model.Run(solution);
  EXPECT_DOUBLE_EQ(solution[0], 1);
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("ConfigVerbose", "0")
                     .WithConfigEntry("NewtonVerbose", "0")
                     .WithConfigEntry("LinearVerbose", "0")
                     .WithConfigEntry("MeshVerbose", "0")
                     .WithConfigEntry("AssembleVerbose", "0")
                     .WithRandomInitialized()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}