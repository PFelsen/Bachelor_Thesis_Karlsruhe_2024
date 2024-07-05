#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "OtherTools.hpp"

TEST(PriorChoice, Adaptive) {
  const int num_iterations = 50;

  std::vector<std::vector<double>> proposal_errors(num_iterations);
  std::string mh_proposals[] = {"RandomWalk", "pCN", "PriorDraw"};

  double precision = 10;
  double variance = 1000;
  double posterior_mean = (1/precision)/((1/precision) + 1/variance);

  for (int i = 0; i < num_iterations; ++i) {
    for (std::string proposal : mh_proposals) {
      SequentialMonteCarlo mcmc{
        SequentialMonteCarloConfig()
          .WithVersion("MCMC")
          .WithNIntermediate(10000)
          .WithNStepsAtOnce(1)
          .WithParticleSetConfig(ParticleSetConfig()
            .WithNParticles(1)
            .WithParticleConfig(ParticleConfig()
              .WithMeasurement({1.0})
              .WithPrecision({precision})
              .WithPDEModelConfig(PDESolverConfig()
                .WithModel("Test")
                .WithProposalConfig(ProposalConfig()
                  .WithProposalType(proposal)
                  .WithRandomField("Simple2DKL")
                  .WithStepLength(0.02)
                  .WithAdaptive()
                  .WithPriorConfig(PriorSamplerConfig()
                    .WithGenerator("Normal")
                    .WithVariance(variance))
                  .WithInitialSampler(PriorSamplerConfig()
                    .WithMean(0.5)
                    .WithVariance(0.5)
                    .WithGenerator("Normal")
                    )))))};
        mcmc.Method();
        std::vector<RVector> markov_chain = mcmc.target_dist.EvaluateParticle(0, 100);
        RVector estimate = Mean(markov_chain);
        proposal_errors[i].push_back( RelativeError1(estimate, RVector(posterior_mean, 1), RVector(1.0, 1)) );
    }
  }
  std::vector<double> proposal_errors_mean = mean_in_column(proposal_errors);
  std::vector<double> proposal_errors_var = variance_in_column(proposal_errors, proposal_errors_mean);

  mout.PrintInfo("SMCHyperparameterData", 1,
    // PrintInfoEntry("Study", "NonAdaptivePriorChoice"),
    PrintInfoEntry("RandomWalkMean", proposal_errors_mean[0]), 
    PrintInfoEntry("RandomWalkVar", proposal_errors_var[0]),
    PrintInfoEntry("pCNMean", proposal_errors_mean[1]),
    PrintInfoEntry("pCNVar", proposal_errors_var[1]),
    PrintInfoEntry("PriorDrawMean", proposal_errors_mean[2]),
    PrintInfoEntry("PriorDrawVar", proposal_errors_var[2]));
}

int main(int argc, char **argv) {
  return MppTest(
    MppTestBuilder(argc, argv).
      WithConfigEntry("ConfigVerbose", "1").
      WithConfigEntry("NewtonVerbose", "0").
      WithConfigEntry("LinearVerbose", "0").
      WithRandomInitialized().
      WithPPM()
  ).RUN_ALL_MPP_TESTS();
}