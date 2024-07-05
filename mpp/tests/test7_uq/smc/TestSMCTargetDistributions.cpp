#include "TestEnvironment.hpp"
#include "TargetDistributions.hpp"

template <typename DistributionType>
void TestInverts(const std::map<std::string, std::vector<double>>& params, double x, double x_test) {
  DistributionType distribution(params);
  double inverted = distribution.CDF(x);
  double double_inverted = distribution.CDFInverse(distribution.CDF(x));
  EXPECT_NEAR(double_inverted, x_test, 1e-6);
}

TEST(ProbabilityDistributions, UniformInverts) {
  std::map<std::string, std::vector<double>> params{{"bounds", {0.0, 1.0}}};
  TestInverts<UniformDistribution>(params, 0.23, 0.23);
  TestInverts<UniformDistribution>(params, -1, 0.0);
  TestInverts<UniformDistribution>(params, 2, 1);
}

TEST(ProbabilityDistributions, NormalInverts) {
  std::map<std::string, std::vector<double>> params{{"mean", {0.0}}, {"standard_deviation", {1.0}}};
  TestInverts<NormalDistribution>(params, 0.23, 0.23);
  TestInverts<NormalDistribution>(params, 4, 4);
}

TEST(ProbabilityDistributions, OrderedMultinomialInverts) {
  std::map<std::string, std::vector<double>> params{{"probabilities", {0.2, 0.1, 0.7}}, {"classes", {0.1, 0.2, 0.3}}};
  TestInverts<OrderedMultinomialDistribution>(params, 0.23, 0.2);
  TestInverts<OrderedMultinomialDistribution>(params, 0.2, 0.2);
  TestInverts<OrderedMultinomialDistribution>(params, -1, 0.1);
  TestInverts<OrderedMultinomialDistribution>(params, 1, 0.3);
}

TEST(ProbabilityDistributions, GetDistribution) {
  std::map<std::string, std::vector<double>> params{{"bounds", {0.0, 1.0}}};
  std::shared_ptr<TargetDistribution> distribution = CreateTargetDistributionShared(
    TargetDistributionConfig()
      .WithParams(params)
      .Type("Uniform"));
  ASSERT_NEAR(distribution->CDF(0.5), 0.5, 1e-6);
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