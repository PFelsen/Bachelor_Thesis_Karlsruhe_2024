#include "TestEnvironment.hpp"
#include "QuantityOfInterest.hpp"

TEST(WeightedSampleDistributionsTest, EvaluateIteration) {
    WeightedSampleDistributions dists;
    dists.AddSample({RVector(1, 1), RVector(2, 1)}, {0.0, 1.0});

    auto [samples, weights] = dists.EvaluateIteration(0);
    
    // Check the samples
    ASSERT_EQ(samples.size(), 2);
    EXPECT_EQ(samples[0], RVector(1, 1));
    EXPECT_EQ(samples[1], RVector(2, 1));
    
    // Check the weights
    ASSERT_EQ(weights.size(), 2);
    EXPECT_DOUBLE_EQ(weights[0], 0.0);
    EXPECT_DOUBLE_EQ(weights[1], 1.0);
}

TEST(WeightedSampleDistributionsTest, EvaluateParticle) {
    WeightedSampleDistributions dists;
    dists.AddSample({RVector(1, 1), RVector(2, 1)});
    dists.AddSample({RVector(3, 1), RVector(4, 1)});
    dists.AddSample({RVector(5, 1), RVector(6, 1)});

    int burnin = 0;
    auto particleSamples = dists.EvaluateParticle(0, burnin);

    // Check the particle samples
    ASSERT_EQ(particleSamples.size(), 3);
    EXPECT_EQ(particleSamples[0], RVector(1, 1));
    EXPECT_EQ(particleSamples[1], RVector(3, 1));
    EXPECT_EQ(particleSamples[2], RVector(5, 1));
}

TEST(WeightedSampleDistributionsTest, EvaluateAll) {
    WeightedSampleDistributions dists;
    dists.AddSample({RVector(1, 1), RVector(2, 1)}, {0.0, 1.0});
    dists.AddSample({RVector(3, 1), RVector(4, 1)}, {0.5, 0.5});
    dists.AddSample({RVector(5, 1), RVector(6, 1)}, {0.5, 0.5});

    int burn_in = 1;
    std::vector<double> iteration_weights = {0.1, 0.9};
    auto [combinedSamples, combinedWeights] = dists.EvaluateAll(burn_in, iteration_weights);

    // Check the combined samples
    ASSERT_EQ(combinedSamples.size(), 4);
    EXPECT_EQ(combinedSamples[0], RVector(3, 1));
    EXPECT_EQ(combinedSamples[1], RVector(4, 1));
    EXPECT_EQ(combinedSamples[2], RVector(5, 1));
    EXPECT_EQ(combinedSamples[3], RVector(6, 1));

    // Check the combined weights
    ASSERT_EQ(combinedWeights.size(), 4);
    EXPECT_DOUBLE_EQ(combinedWeights[0], 0.5 * 0.1);
    EXPECT_DOUBLE_EQ(combinedWeights[1], 0.5 * 0.1);
    EXPECT_DOUBLE_EQ(combinedWeights[2], 0.5 * 0.9);
    EXPECT_DOUBLE_EQ(combinedWeights[3], 0.5 * 0.9);
}

TEST(WeightedSampleDistributionsTest, CalculateWeightedMean) {
    // Test for unweighted mean
    std::vector<RVector> unweightedSamples = {RVector{1, 2}, RVector{3, 4}, RVector{5, 6}};
    auto unweightedMean = Mean(unweightedSamples);
    EXPECT_DOUBLE_EQ(unweightedMean[0], 3.0);
    EXPECT_DOUBLE_EQ(unweightedMean[1], 4.0);

    // Test for weighted mean
    std::vector<RVector> weightedSamples = {RVector{1, 2}, RVector{3, 4}, RVector{5, 6}};
    std::vector<double> weights = {0.5, 0.3, 0.2};
    auto weightedMean = Mean(weightedSamples, weights);
    EXPECT_DOUBLE_EQ(weightedMean[0], (1*0.5 + 3*0.3 + 5*0.2));
    EXPECT_DOUBLE_EQ(weightedMean[1], (2*0.5 + 4*0.3 + 6*0.2));
}

TEST(TotalError1Test, ComputesCorrectly) {
    RVector estimate(std::vector<double>{0, 1});
    RVector true_parameter(std::vector<double>{-1, 1});
    RVector sqrtEigenvalues(std::vector<double>{1, 0.5});

    double total_error = TotalError1(estimate, true_parameter, sqrtEigenvalues);
    EXPECT_DOUBLE_EQ(total_error, 1.0);
}

TEST(RelativeError1Test, ComputesCorrectly) {
    RVector estimate(std::vector<double>{0, 1});
    RVector true_parameter(std::vector<double>{-1, 1});
    RVector sqrtEigenvalues(std::vector<double>{1, 0.5});

    double relative_error = RelativeError1(estimate, true_parameter, sqrtEigenvalues);
    EXPECT_DOUBLE_EQ(relative_error, 2.0 / 3.0);
}

TEST(GoodnessOfFit, TwoSampleKLDistance) {
  std::vector<double> sample1, sample2, weights1, weights2;
  sample1 = {0.1, 0.2};
  sample2 = {0.1, 0.15};
  weights2 = {0.75, 0.25};
  EXPECT_NEAR(TwoSampleKLDistance(sample1, sample2, weights2), 0.5, 1e-6);
  sample2 = {0.05, 0.15};
  EXPECT_NEAR(TwoSampleKLDistance(sample1, sample2, weights2), 0.75, 1e-6);
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