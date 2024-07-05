#include "TestEnvironment.hpp"
#include "Metropolis.hpp"
#include "PDEModel.hpp"
#include "KolmogorovSmirnov.hpp"

TEST(Metropolis, TestMeasurement) {
  std::shared_ptr<Measurement> measurement = std::make_shared<Measurement>(
    RVector({1.0}), RMatrix({{1.0}})
  );
  EXPECT_DOUBLE_EQ(measurement->LogLikelihood(RVector(0.5, 1)), 0.125);
}

TEST(Metropolis, ReplicatedStep) {
  RVector proposal(1), solution(1);
  ProposalConfig prop_conf = ProposalConfig()
    .WithProposalType("RandomWalk")
    .WithStepLength(0.5)
    .WithRandomField("Simple2DKL")
    .WithPriorConfig(PriorSamplerConfig().WithGenerator("Normal").WithVariance({1}));
  std::unique_ptr<ProposalGenerator> proposal_kernel =  CreateProposalGeneratorUnique(prop_conf, proposal);
  proposal_kernel->GetState()[0] = 0.5;
  proposal_kernel->GetProposal()[0] = 1;

  TestModel pde_model(PDESolverConfig().WithProposalConfig(prop_conf).WithProblem("Test1D"));
  pde_model.GetProblem().GetInput().GetRandomField().GetParameterVector() = proposal_kernel->GetProposal();
  pde_model.Run(solution);
  EXPECT_DOUBLE_EQ(solution[0], 1);

  Measurement measurement({0}, {{1}});
  double old_potential = measurement.LogLikelihood(proposal_kernel->GetState());
  EXPECT_DOUBLE_EQ(old_potential, 0.125);
  double new_potential = measurement.LogLikelihood(solution);
  EXPECT_DOUBLE_EQ(new_potential, 0.5);

  double acceptance_probability = proposal_kernel->AcceptanceProbability(old_potential, new_potential);
  EXPECT_DOUBLE_EQ(acceptance_probability, exp(-0.75));

  double random_draw = 0.0;
  if(acceptance_probability > random_draw) {
    proposal_kernel->Accept();
  }
  EXPECT_DOUBLE_EQ(proposal_kernel->GetState()[0], 1);
}

TEST(Metropolis, ManualMetropolisMarkovChain) {
  RVector solution(1);
  ProposalConfig prop_conf = ProposalConfig()
    .WithProposalType("RandomWalk")
    .WithStepLength(2)
    .WithRandomField("Simple2DKL")
    .WithPriorConfig(PriorSamplerConfig().WithGenerator("Normal").WithVariance({1}));
  TestModel pde_model(PDESolverConfig().WithProposalConfig(prop_conf).WithProblem("Test1D"));
  std::unique_ptr<ProposalGenerator> proposal_kernel = 
    CreateProposalGeneratorUnique(prop_conf, pde_model.GetProblem().GetInput().GetRandomField().GetParameterVector());
  Measurement measurement({1}, {{1}});
  double old_potential = measurement.LogLikelihood(proposal_kernel->GetState());
  double new_potential, acceptance_probability, random_draw;

  std::vector<double> acceptance_probabilities;
  std::vector<double> markov_chain;
  for(int j = 0; j < 1000; j++) {
    proposal_kernel->Propose();
    pde_model.Run(solution);
    new_potential = measurement.LogLikelihood(solution);

    acceptance_probability = proposal_kernel->AcceptanceProbability(old_potential, new_potential);

    random_draw = Random::Uniform(0, 0, 1);
    if(acceptance_probability > random_draw) {
      proposal_kernel->Accept();
      old_potential = new_potential;
    }
    acceptance_probabilities.push_back(acceptance_probability);
    if(j > 200 && j%8 ==0) markov_chain.push_back(proposal_kernel->GetState()[0]);
  }


  KolmogorovSmirnovTest tester(
    TargetDistributionConfig().
      Type("Normal").
      WithParams({{"mean", {0.5}}, {"variance", {0.5}}}));
  EXPECT_TRUE(tester.Run(markov_chain, 0.05));
}

void TestMetropolis(const std::string& proposalType) {
    ProposalConfig proposal_config = ProposalConfig()
        .WithStepLength(2)
        .WithProposalType(proposalType)
        .WithRandomField("Simple2DKL")
        .WithPriorConfig(PriorSamplerConfig()
            .WithGenerator("Normal")
            .WithVariance(1.0));

    std::shared_ptr<EllipticPDEModel> model = std::make_shared<TestModel>(
        PDESolverConfig().WithProposalConfig(proposal_config).WithProblem("Test1D"));

    std::shared_ptr<Measurement> measurement = std::make_shared<Measurement>(
        RVector({1.0}), RMatrix({{1.0}})
    );

    MetropolisKernel kernel(model, measurement, proposal_config);

    std::vector<double> markov_chain;

    for (int j = 0; j < 1000; j++) {
        kernel.Step();
        if (j > 50 && j % 5 == 0) {
            markov_chain.push_back(kernel.GetProposalKernel().GetState()[0]);
        }
    }

    KolmogorovSmirnovTest tester(
        TargetDistributionConfig()
            .Type("Normal")
            .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
    );

    EXPECT_TRUE(tester.Run(markov_chain, 0.05));
}

TEST(Metropolis, RandomWalkMetropolis) { TestMetropolis("RandomWalk"); }
TEST(Metropolis, pCNMetropolis) { TestMetropolis("pCN"); }
TEST(Metropolis, PriorDrawMetropolis) { TestMetropolis("PriorDraw"); }

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