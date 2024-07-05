#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "OtherTools.hpp"

TEST(SMC, ConvergenceInM) {
  const int num_iterations = 50;

  std::vector<std::vector<double>> errors(num_iterations);
  std::vector<int> increasingM = {5, 10, 20, 40, 70, 100, 200, 400, 700, 1000, 2000, 4000, 7000};

  double precision = 10;
  double variance = 1000;
  double posterior_mean = (1/precision)/((1/precision) + 1/variance);

  for (int i = 0; i < num_iterations; ++i) {
    for (int M : increasingM) {
      SequentialMonteCarlo smc{
        SequentialMonteCarloConfig()
          .WithVersion("SMC")
          .WithNIntermediate(20)
          .WithNStepsAtOnce(1)
          .WithRelativeESSMin(0.8)
          .WithParticleSetConfig(ParticleSetConfig()
            .WithNParticles(M)
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
    // PrintInfoEntry("Study", "ConvergenceInM"),
    PrintInfoEntry("SMCErrorMean", errors_mean), 
    PrintInfoEntry("SMCErrorVar", errors_var),
    PrintInfoEntry("MParticles", increasingM));
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