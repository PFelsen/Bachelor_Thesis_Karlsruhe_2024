#include "TestEnvironment.hpp"
#include "KolmogorovSmirnov.hpp"
#include "ChiSquared.hpp"
#include "ShapiroWilk.hpp"
#include "Binning.hpp"

// Many of these tests derived by comparison with implementations from python scipy.stats

TEST(GoodnessOfFit, KolmovorovSmirnovPValue) {
  KolmogorovSmirnovTest kstest{TargetDistributionConfig()};
  EXPECT_NEAR(kstest.p_value(0.9, 1), 0.2, 1e-3);
  EXPECT_NEAR(kstest.p_value(0.270, 35), 0.01, 1e-3);
  EXPECT_NEAR(kstest.p_value(0.11139799127367866, 100), 0.1548, 1e-3);
}

TEST(GoodnessOfFit, Chi2PValue) {
  ChiSquaredTest cstest{TargetDistributionConfig()};
  EXPECT_NEAR(cstest.p_value(23.685, 14), 0.05, 1e-3);
  EXPECT_NEAR(cstest.p_value(1, 6), 0.9857, 1e-3);
  EXPECT_NEAR(cstest.p_value(20, 6), 0.0028, 1e-3);
  EXPECT_NEAR(cstest.p_value(170, 130), 0.0106, 1e-3);
  EXPECT_NEAR(cstest.p_value(190, 130), 0.0, 1e-3);
}

TEST(GoodnessOfFit, KolmovorovSmirnovStatistic) {
  std::vector<double> sample;
  std::map<std::string, std::vector<double>> params{{"bounds", {0.0, 1.0}}};
  KolmogorovSmirnovTest kstest{
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Uniform")};
  sample = {0.5};
  EXPECT_NEAR(kstest.test_statistic(sample), 0.5, 1e-6);
  sample = {-10.0};
  EXPECT_NEAR(kstest.test_statistic(sample), 1.0, 1e-6);
  sample = {0.1, 0.3, 0.5, 0.7, 0.9};
  EXPECT_NEAR(kstest.test_statistic(sample), 0.1, 1e-6);
}

TEST(GoodnessOfFit, Chi2Statistic) {
  ChiSquaredTest cstest{TargetDistributionConfig()};
  EXPECT_NEAR(cstest.chi2_statistic(
    {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, 
    {5, 5, 5, 5, 5, 5, 5, 5, 5, 5}), 0.0, 1e-6);
  EXPECT_NEAR(cstest.chi2_statistic(
    {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, 
    {5, 5, 5, 5, 5, 5, 5, 5, 5, 10}), 4.09090909, 1e-6);
}

TEST(GoodnessOfFit, Chi2InverseCDFBinning) {
  std::map<std::string, std::vector<double>> params{{"bounds", {0.0, 1.0}}};
  InverseCDFBinning binning{
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Uniform")};
  int sample_size = 50;
  std::vector<double> solution = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  
  std::vector<double> delimiters = binning.binning_delimiters(sample_size);
  std::vector<double> probabilities = binning.probabilities(delimiters);

  for (size_t i = 0; i < solution.size(); ++i) {
    EXPECT_NEAR(solution[i], delimiters[i], 1e-6)
      << "Element at index " << i << " differs significantly.";
    EXPECT_NEAR(probabilities[i], 0.1, 1e-6);
  }
}

TEST(GoodnessOfFit, Chi2EmpiricalBinning) {
  std::map<std::string, std::vector<double>> params{{"bounds", {0.0, 1.0}}};
  EmpiricalBinning binning{
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Uniform")};
  int sample_size = 70;

  std::vector<double> delimiters = binning.binning_delimiters(sample_size);
  std::vector<double> probabilities = binning.probabilities(delimiters);

  EXPECT_GE(delimiters[0], 5.0/sample_size);
  EXPECT_GE(1 - delimiters.back(), 5.0/sample_size);
  for (int i = 0; i < delimiters.size() - 1; i++) {
    EXPECT_GE(delimiters[i+1] - delimiters[i], 5.0/sample_size);
    EXPECT_GE(delimiters[i+1] - delimiters[i], 5.0/sample_size);
    EXPECT_GE(probabilities[i], 5.0/sample_size);
  }
}

TEST(GoodnessOfFit, Chi2Histogram) {
  ChiSquaredTest cstest{TargetDistributionConfig()};
  std::vector<int> solution = {0,3,0};
  EXPECT_EQ(cstest.count_bin_samples({0.8, 0.6, 1}, {0.5, 1}), solution);
}

TEST(GoodnessOfFit, KolmogorovSmirnovWeighted) {
  std::vector<double> sample, weights;
  std::map<std::string, std::vector<double>> params{{"bounds", {0.0, 1.0}}};
  KolmogorovSmirnovTest kstest{
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Uniform")};
  sample = {0.3, 0.6};
  weights = {0.5, 0.5};
  EXPECT_NEAR(kstest.test_statistic(sample), 0.4, 1e-6);
  sample = {0.3, 0.5};
  weights = {0.4, 0.6};
  EXPECT_NEAR(kstest.test_statistic(sample), 0.5, 1e-6);
}

TEST(GoodnessOfFit, KolmogorovSmirnovTwoSample) {
  std::vector<double> sample1, sample2, weights1, weights2;
  weights1 = {0.2, 0.1, 0.7};
  sample1 = {0.1, 0.2, 0.3};
  std::map<std::string, std::vector<double>> params{{"probabilities", weights1}, {"classes", sample1}};
  KolmogorovSmirnovTest kstest{
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Multinomial")};
  sample2 = {0.05, 0.25, 0.5};
  weights2 = {0.35, 0.25, 0.4};
  EXPECT_NEAR(kstest.test_statistic(sample2, weights2), 0.4, 1e-6);
}

/*
TEST(TestGoodnessOfFitTesters, SmapiroWilkPValue) {
  ShapiroWilkTest shapiro{TargetDistributionConfig()};
  EXPECT_NEAR(shapiro.p_value(100, 0.9884352684020996), 0.5408, 1e-3);
}

TEST(TestGoodnessOfFitTesters, ShapiroWilkStatistic) {
  std::vector<double> sample;
  std::map<std::string, std::vector<double>> params{{"mean", {0.0}}, {"standard_deviation", {1.0}}};
  KolmogorovSmirnovTest kstest{
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Normal")};
  sample = {-1.0, 0.0, 1.0};
  EXPECT_NEAR(kstest.test_statistic(sample), 1.0, 1e-3);
}
*/

TEST(SortSampleAndWeightsTest, SortAndRearrange) {
    // Test case data
    std::vector<double> sample = {3.0, 1.0, 2.0};
    std::vector<double> weights = {0.3, 0.1, 0.2};

    // Expected results after sorting and rearranging
    std::vector<double> expectedSample = {1.0, 2.0, 3.0};
    std::vector<double> expectedWeights = {0.1, 0.2, 0.3};

    // Call the function to sort and rearrange
    SortSampleAndWeights(sample, weights);

    // Check if the sample and weights are as expected
    ASSERT_EQ(sample, expectedSample);
    ASSERT_EQ(weights, expectedWeights);
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