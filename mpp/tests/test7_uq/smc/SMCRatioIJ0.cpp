#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "OtherTools.hpp"

TEST(SMC, BestIJRatio) {
  const int num_iterations = 50;

  std::vector<std::vector<double>> errors(num_iterations);
  std::vector<int> num_intermediate = {1, 2, 4, 5, 10, 20, 25, 50, 100};

  double precision = 10;
  double variance = 1000;
  double posterior_mean = (1/precision)/((1/precision) + 1/variance);
  int J;

  for (int i = 0; i < num_iterations; ++i) {
    for (int I : num_intermediate) {
      J = 100/I;
      SequentialMonteCarlo smc{
        SequentialMonteCarloConfig()
          .WithVersion("SMC")
          .WithNIntermediate(I)
          .WithNStepsAtOnce(J)
          .WithRelativeESSMin(0.5)
          .WithParticleSetConfig(ParticleSetConfig()
            .WithNParticles(100)
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
    PrintInfoEntry("IJRatioErrorMean", errors_mean),
    PrintInfoEntry("IJRatioErrorVar", errors_var),
    PrintInfoEntry("NIntermediateMeasures", num_intermediate));
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