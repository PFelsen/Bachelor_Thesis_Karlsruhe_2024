#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "KolmogorovSmirnov.hpp"

TEST(TestSMC, AsMCMC) {
  SequentialMonteCarlo mcmc{
    SequentialMonteCarloConfig()
      .WithVersion("MCMC")
      .WithNIntermediate(1000)
      .WithNStepsAtOnce(4)
      .WithParticleSetConfig(ParticleSetConfig()
        .WithResamplingScheme("Multinomial")
        .WithNParticles(1)
        .WithParticleConfig(ParticleConfig()
          .WithMeasurement({1.0})
          .WithPrecision({1.0})
          .WithPDEModelConfig(PDESolverConfig()
            .WithModel("Test")
            .WithProposalConfig(ProposalConfig()
              .WithProposalType("RandomWalk")
              .WithRandomField("Simple2DKL")
              .WithStepLength(2)
              .WithPriorConfig(PriorSamplerConfig()
                .WithGenerator("Normal")
                .WithVariance(1.0)))
              )))};
  mcmc.Method();
  int burn_in = 10;
  std::vector<double> solution = ExtractCoordinateData(mcmc.target_dist.EvaluateParticle(0, burn_in), 0);

  KolmogorovSmirnovTest tester{
    TargetDistributionConfig()
      .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
      .Type("Normal")};
  EXPECT_TRUE(tester.Run(solution, 0.05));
}

TEST(TestSMC, AsImportanceSampling) {
  SequentialMonteCarlo smc{
    SequentialMonteCarloConfig()
      .WithVersion("IS")
      .WithParticleSetConfig(ParticleSetConfig()
        .WithResamplingScheme("Multinomial")
        .WithNParticles(1000)
        .WithParticleConfig(ParticleConfig()
          .WithMeasurement({1.0})
          .WithPrecision({1.0})
          .WithPDEModelConfig(PDESolverConfig()
            .WithModel("Test")
            .WithProposalConfig(ProposalConfig()
              .WithProposalType("RandomWalk")
              .WithRandomField("Simple2DKL")
              .WithStepLength(2)
              .WithPriorConfig(PriorSamplerConfig()
                .WithGenerator("Normal")
                .WithVariance(1.0)))
              )))};
  smc.Method();
  std::tuple<std::vector<RVector>, std::vector<double>> data = smc.target_dist.EvaluateIteration(1);
  std::vector<double> sample = ExtractCoordinateData(std::get<0>(data), 0);
  std::vector<double> weights = std::get<1>(data);

  KolmogorovSmirnovTest tester{
    TargetDistributionConfig()
      .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
      .Type("Normal")};
  EXPECT_TRUE(tester.Run(sample, 0.001, weights));
}

TEST(TestSMC, SMCNoMarkovStep) {
  int n_intermediate = 10;
  SequentialMonteCarlo smc{
    SequentialMonteCarloConfig()
      .WithNIntermediate(n_intermediate)
      .WithNStepsAtOnce(0)
      .WithRelativeESSMin(0)
      .WithVersion("SMC")
      .WithParticleSetConfig(ParticleSetConfig()
        .WithResamplingScheme("Multinomial")
        .WithNParticles(1000)
        .WithParticleConfig(ParticleConfig()
          .WithMeasurement({1.0})
          .WithPrecision({1.0})
          .WithPDEModelConfig(PDESolverConfig()
            .WithModel("Test")
            .WithProposalConfig(ProposalConfig()
              .WithRandomField("Simple2DKL")
              .WithPriorConfig(PriorSamplerConfig()
                .WithGenerator("Normal")
                .WithVariance(1.0)))
              )))};
  smc.Method();
  std::tuple<std::vector<RVector>, std::vector<double>> data = smc.target_dist.EvaluateIteration(n_intermediate+1);
  std::vector<double> sample = ExtractCoordinateData(std::get<0>(data), 0);
  std::vector<double> weights = std::get<1>(data);

  KolmogorovSmirnovTest tester{
    TargetDistributionConfig()
      .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
      .Type("Normal")};
  EXPECT_TRUE(tester.Run(sample, 0.05, weights));
}

TEST(TestSMC, SMCOnlyMarkovSteps) {
  SequentialMonteCarlo smc{
    SequentialMonteCarloConfig()
      .WithNIntermediate(0)
      .WithNStepsAtOnce(10)
      .WithRelativeESSMin(0)
      .WithVersion("SMC")
      .WithParticleSetConfig(ParticleSetConfig()
        .WithResamplingScheme("Multinomial")
        .WithNParticles(1000)
        .WithParticleConfig(ParticleConfig()
          .WithMeasurement({1.0})
          .WithPrecision({1.0})
          .WithPDEModelConfig(PDESolverConfig()
            .WithModel("Test")
            .WithProposalConfig(ProposalConfig()
              .WithRandomField("Simple2DKL")
              .WithPriorConfig(PriorSamplerConfig()
                .WithGenerator("Normal")
                .WithVariance(1.0)))
              )))};
  smc.Method();
  std::tuple<std::vector<RVector>, std::vector<double>> data = smc.target_dist.EvaluateIteration(1);
  std::vector<double> sample = ExtractCoordinateData(std::get<0>(data), 0);
  std::vector<double> weights = std::get<1>(data);

  KolmogorovSmirnovTest tester{
    TargetDistributionConfig()
      .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
      .Type("Normal")};
  EXPECT_TRUE(tester.Run(sample, 0.05, weights));
}

TEST(TestSMC, WithoutResampling) {
  int n_intermediate = 4;
  SequentialMonteCarlo smc{
    SequentialMonteCarloConfig()
      .WithNIntermediate(n_intermediate)
      .WithNStepsAtOnce(4)
      .WithRelativeESSMin(0)
      .WithVersion("SMC")
      .WithParticleSetConfig(ParticleSetConfig()
        .WithResamplingScheme("Multinomial")
        .WithNParticles(1000)
        .WithParticleConfig(ParticleConfig()
          .WithMeasurement({1.0})
          .WithPrecision({1.0})
          .WithPDEModelConfig(PDESolverConfig()
            .WithModel("Test")
            .WithProposalConfig(ProposalConfig()
              .WithRandomField("Simple2DKL")
              .WithPriorConfig(PriorSamplerConfig()
                .WithGenerator("Normal")
                .WithVariance(1.0)))
              )))};
  smc.Method();
  std::tuple<std::vector<RVector>, std::vector<double>> data = smc.target_dist.EvaluateIteration(n_intermediate + 1);
  std::vector<double> sample = ExtractCoordinateData(std::get<0>(data), 0);
  std::vector<double> weights = std::get<1>(data);

  KolmogorovSmirnovTest tester{
    TargetDistributionConfig()
      .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
      .Type("Normal")};
//  TODO: FIX ME
//  EXPECT_TRUE(tester.Run(sample, 0.05, weights));
}

TEST(TestSMC, TryToResample) {
  int n_intermediate = 10;
  SequentialMonteCarlo smc{
    SequentialMonteCarloConfig()
      .WithNIntermediate(n_intermediate)
      .WithRelativeESSMin(1)
      .WithNStepsAtOnce(5)
      .WithVersion("SMC")
      .WithParticleSetConfig(ParticleSetConfig()
        .WithResamplingScheme("Multinomial")
        .WithNParticles(100)
        .WithParticleConfig(ParticleConfig()
          .WithMeasurement({1.0})
          .WithPrecision({1.0})
          .WithPDEModelConfig(PDESolverConfig()
            .WithModel("Test")
            .WithProposalConfig(ProposalConfig()
              .WithRandomField("Simple2DKL")
              .WithPriorConfig(PriorSamplerConfig()
                .WithGenerator("Normal")
                .WithVariance(1.0)))
              )))};
  smc.Method();
  std::tuple<std::vector<RVector>, std::vector<double>> data = smc.target_dist.EvaluateIteration(n_intermediate+1);
  std::vector<double> sample = ExtractCoordinateData(std::get<0>(data), 0);
  std::vector<double> weights = std::get<1>(data);

  KolmogorovSmirnovTest tester{
    TargetDistributionConfig()
      .WithParams({{"mean", {0.5}}, {"variance", {0.5}}})
      .Type("Normal")};
  EXPECT_TRUE(tester.Run(sample, 0.01, weights));
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