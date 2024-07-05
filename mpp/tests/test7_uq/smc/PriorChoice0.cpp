#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "OtherTools.hpp"

TEST(PriorChoice, Adaptive) {
  const int num_iterations = 50;
  double precision = 10;

  std::vector<std::vector<double>> errors_over_variance_scaling(num_iterations);
  std::vector<double> variance_scale;

  for (int i = 0; i < num_iterations; ++i) {
    for (int j = 0; j < 20; j++) {
      if (i == 0) variance_scale.push_back(std::pow(2, j));
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
                  .WithProposalType("RandomWalk")
                  .WithRandomField("Simple2DKL")
                  .WithStepLength(0.02)
                  .WithAdaptive()
                  .WithPriorConfig(PriorSamplerConfig()
                    .WithGenerator("Normal")
                    .WithVariance(std::pow(2, j)))
                  .WithInitialSampler(PriorSamplerConfig()
                    .WithMean((1/precision)/((1/precision) + 1/std::pow(2, j)))
                    .WithVariance(1/(1/precision + 1/std::pow(2, j)))
                    .WithGenerator("Normal")
                    )))))};
        mcmc.Method();
        std::vector<RVector> markov_chain = mcmc.target_dist.EvaluateParticle(0, 0);
        RVector estimate = Mean(markov_chain);
        errors_over_variance_scaling[i].push_back( RelativeError1(estimate, RVector(1.0, 1), RVector(1.0, 1)) );
    }
  }
  std::vector<double> errors_over_variance_mean = mean_in_column(errors_over_variance_scaling);
  std::vector<double> errors_over_variance_var = variance_in_column(errors_over_variance_scaling, errors_over_variance_mean);

  mout.PrintInfo("SMCHyperparameterData", 1, 
    // PrintInfoEntry("Study", "PriorChoiceNotAdaptive"),
    PrintInfoEntry("MeanErrorsOverPriorVarianceScale", errors_over_variance_mean), 
    PrintInfoEntry("VarErrorsOverPriorVarianceScale", errors_over_variance_var), 
    PrintInfoEntry("PriorVarianceScale", variance_scale));
}

TEST(PriorChoice, NotAdaptive) {
  const int num_iterations = 50;
  double precision = 10;

  std::vector<std::vector<double>> errors_over_variance_scaling(num_iterations);
  std::vector<double> variance_scale;

  for (int i = 0; i < num_iterations; ++i) {
    for (int j = 0; j < 20; j++) {
      if (i == 0) variance_scale.push_back(std::pow(2, j));
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
                  .WithProposalType("RandomWalk")
                  .WithRandomField("Simple2DKL")
                  .WithStepLength(0.02)
                  .WithPriorConfig(PriorSamplerConfig()
                    .WithGenerator("Normal")
                    .WithVariance(std::pow(2, j)))
                  .WithInitialSampler(PriorSamplerConfig()
                    .WithMean((1/precision)/((1/precision) + 1/std::pow(2, j)))
                    .WithVariance(1/(1/precision + 1/std::pow(2, j)))
                    .WithGenerator("Normal")
                    )))))};
        mcmc.Method();
        std::vector<RVector> markov_chain = mcmc.target_dist.EvaluateParticle(0, 0);
        RVector estimate = Mean(markov_chain);
        errors_over_variance_scaling[i].push_back( RelativeError1(estimate, RVector(1.0, 1), RVector(1.0, 1)) );
    }
  }
  std::vector<double> errors_over_variance_mean = mean_in_column(errors_over_variance_scaling);
  std::vector<double> errors_over_variance_var = variance_in_column(errors_over_variance_scaling, errors_over_variance_mean);

  mout.PrintInfo("BurnInAndMarkovChainLength", 1, 
    // PrintInfoEntry("Study", "PriorChoiceNotAdaptive"),
    PrintInfoEntry("MeanErrorsOverPriorVarianceScaleNotAdaptive", errors_over_variance_mean), 
    PrintInfoEntry("VarErrorsOverPriorVarianceScaleNotAdaptive", errors_over_variance_var), 
    PrintInfoEntry("PriorVarianceScaleNotAdaptive", variance_scale));
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