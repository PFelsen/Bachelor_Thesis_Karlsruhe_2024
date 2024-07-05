#include "ChiSquared.hpp"
#include "KolmogorovSmirnov.hpp"
#include "OtherTools.hpp"
#include "Random.hpp"
#include "Resampling.hpp"
#include "TestEnvironment.hpp"

TEST(Resampler, ParticleCopying) {
  std::vector<double> weights = {1, 0};
  KolmogorovSmirnovTest tester{
      TargetDistributionConfig().WithParams({{"probabilities", weights}}).Type("Multinomial")};
  std::unique_ptr<Resampling> resampler = CreateResamplingUnique("Multinomial");
  std::shared_ptr<Particle> particle1 =
      std::make_shared<Particle>(ParticleConfig().WithPDEModelConfig(
          PDESolverConfig().WithProblem("Test1D").WithModel("Lagrange")));
  std::shared_ptr<Particle> particle2 =
      std::make_shared<Particle>(ParticleConfig().WithPDEModelConfig(
          PDESolverConfig().WithProblem("Test1D").WithModel("Lagrange")));
  particle1->GetKernel().GetProposalKernel().SetState({10});
  std::vector<std::shared_ptr<Particle>> particlevec = {particle1, particle2};
  EXPECT_EQ(particlevec[0]->GetKernel().GetProposalKernel().GetState()[0], 10);
  std::vector<int> indices = resampler->getIndices(weights);
  EXPECT_EQ(indices[0], 0);
  EXPECT_EQ(indices[1], 0);
  resampler->resample(weights, particlevec);
  EXPECT_EQ(particle2->GetKernel().GetProposalKernel().GetState()[0], 10);
}

TEST(Resampler, MultinomialEqualProbabilities) {
  std::vector<double> weights(40, 0.025);
  KolmogorovSmirnovTest tester{
      TargetDistributionConfig().WithParams({{"probabilities", weights}}).Type("Multinomial")};
  std::unique_ptr<Resampling> resampler = CreateResamplingUnique("Multinomial");
  std::vector<double> sample;
  for (int i = 0; i < 25; i++) {
    for (int index : resampler->getIndices(weights)) {
      sample.push_back(index);
    }
  }
  EXPECT_TRUE(tester.Run(sample, 0.05));
}

TEST(Resampler, MultinomialUnequalProbabilities) {
  RVector uniform_draw = Random::Uniform(0, 40, 0, 1);
  std::vector<double> weights = uniform_draw.Data();
  normalizeWeights(weights);
  KolmogorovSmirnovTest tester{
      TargetDistributionConfig().WithParams({{"probabilities", weights}}).Type("Multinomial")};
  std::unique_ptr<Resampling> resampler = CreateResamplingUnique("Multinomial");
  std::vector<double> sample;
  for (int i = 0; i < 25; i++) {
    for (int index : resampler->getIndices(weights)) {
      sample.push_back(index);
    }
  }
  EXPECT_TRUE(tester.Run(sample, 0.05));
}

TEST(Resampler, SystematicEqualProbabilities) {
  std::vector<double> weights(40, 0.025);
  KolmogorovSmirnovTest tester{
      TargetDistributionConfig().WithParams({{"probabilities", weights}}).Type("Multinomial")};
  std::unique_ptr<Resampling> resampler = CreateResamplingUnique("Systematic");
  std::vector<double> sample;
  for (int i = 0; i < 25; i++) {
    for (int index : resampler->getIndices(weights)) {
      sample.push_back(index);
    }
  }
  EXPECT_TRUE(tester.Run(sample, 0.05));
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