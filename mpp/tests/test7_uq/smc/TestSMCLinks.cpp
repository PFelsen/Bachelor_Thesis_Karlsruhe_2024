#include "TestEnvironment.hpp"
#include "PDEModel.hpp"
#include "Metropolis.hpp"

void RFParameterTest(const std::string &problemName, const std::string &randomFieldName) {
  EllipticPDEModel model(PDESolverConfig()
                             .WithProblem(problemName)
                             .WithProposalConfig(
                                 ProposalConfig().WithRandomField(randomFieldName)));
  RVector &ref1 = model.GetProblem().GetInput().GetRandomField().GetParameterVector();
  RVector &ref2 = model.GetProblem().GetInput().GetRandomField().GetParameterVector();
  ref1[0] = 10;
  EXPECT_EQ(ref2[0], 10);
}

TEST(LinkingSMC, RandomFieldParameterToPDEModel0) { RFParameterTest("Test2D", "Simple2DKL"); }

TEST(LinkingSMC, RandomFieldParameterToPDEModel1) { RFParameterTest("Test1D", "Exp2DKL"); }

TEST(LinkingSMC, RandomFieldParameterToPDEModel2) { RFParameterTest("RainLaplace2D", "Long2DKL"); }

TEST(LinkingSMC, RandomFieldParameterToProposalGenerator) {
  ProposalConfig proposal_config;
  std::shared_ptr<EllipticPDEModel> model = std::make_shared<EllipticPDEModel>(PDESolverConfig().
      WithProposalConfig(ProposalConfig().
      WithRandomField("Exp2DKL")));
  std::shared_ptr<Measurement> measurement = std::make_shared<Measurement>(
      RVector({1.0}), RMatrix({{1.0}})
  );
  MetropolisKernel metropolis(model, measurement, proposal_config);
  model->GetProblem().GetInput().GetRandomField().GetParameterVector()[0] = 10;
  EXPECT_EQ(metropolis.GetProposalKernel().GetProposal()[0], 10);
  metropolis.GetProposalKernel().Propose();
  EXPECT_EQ(metropolis.GetProposalKernel().GetProposal()[0],
            model->GetProblem().GetInput().GetRandomField().GetParameterVector()[0]);
}

int main(int argc, char **argv) {
  return MppTest(
      MppTestBuilder(argc, argv).
          WithConfigEntry("ConfigVerbose", "0").
          WithConfigEntry("NewtonVerbose", "0").
          WithConfigEntry("LinearVerbose", "0").
          WithConfigEntry("MeshVerbose", "0").
          WithConfigEntry("AssembleVerbose", "0").
          WithRandomInitialized().
          WithPPM()
  ).RUN_ALL_MPP_TESTS();
}