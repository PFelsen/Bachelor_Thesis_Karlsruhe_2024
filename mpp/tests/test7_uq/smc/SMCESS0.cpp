#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "OtherTools.hpp"

TEST(SMC, BestESS) {
  const int num_iterations = 50;

  std::vector<std::vector<double>> errors(num_iterations);
  std::vector<double> ess_min;
  for (double ess = 0; ess < 1; ess+=0.1) ess_min.push_back(ess);

  double precision = 10;
  double variance = 1000;
  double posterior_mean = (1/precision)/((1/precision) + 1/variance);

  for (int i = 0; i < num_iterations; ++i) {
    for (double ess : ess_min) {
      SequentialMonteCarlo smc{
        SequentialMonteCarloConfig()
          .WithVersion("SMC")
          .WithNIntermediate(50)
          .WithNStepsAtOnce(1)
          .WithRelativeESSMin(ess)
          .WithParticleSetConfig(ParticleSetConfig()
            .WithNParticles(1000)
            .WithParticleConfig(ParticleConfig()
              .WithMeasurement({1.0})
              .WithPrecision({precision})
              .WithPDEModelConfig(PDESolverConfig()
                .WithModel("Test")
                .WithProposalConfig(ProposalConfig()
                  .WithProposalType("RandomWalk")
                  .WithRandomField("Simple2DKL")
                  .WithPriorConfig(PriorSamplerConfig()
                    .WithGenerator("Normal")
                    .WithVariance(variance))
                  .WithInitialSampler(PriorSamplerConfig()
                    .WithMean(0.5)
                    .WithVariance(0.5)
                    .WithGenerator("Normal")
                    )))))};
        smc.Method();
        std::tuple<std::vector<RVector>, std::vector<double>> data = smc.target_dist.EvaluateIteration(-1);
        RVector estimate = Mean(std::get<0>(data), std::get<1>(data));
        errors[i].push_back( RelativeError1(estimate, RVector(posterior_mean, 1), RVector(1.0, 1)) );
    }
  }
  std::vector<double> errors_mean = mean_in_column(errors);
  std::vector<double> errors_var = variance_in_column(errors, errors_mean);

  mout.PrintInfo("SMCHyperparameterData", 1, 
    // PrintInfoEntry("Study", "BestESS"),
    PrintInfoEntry("ESSErrorMean", errors_mean), 
    PrintInfoEntry("ESSErrorVar", errors_var),
    PrintInfoEntry("ESSMin", ess_min));
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