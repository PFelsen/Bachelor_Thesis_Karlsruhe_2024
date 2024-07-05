#include "TestRandom.hpp"

TEST_UNIFORM(TestUniformDistribution021)

TEST_UNIFORM(TestUniformDistribution021WithSplit)

TEST_UNIFORM(TestUniformDistribution021WithDoubleSplit)

TEST_UNIFORM(TestUniformDistribution021WithFullSplit)

TEST_UNIFORM(TestUniformDistributionNeg12Pos1)

TEST_UNIFORM(TestUniformDistributionNeg12Pos1WithSplit)

TEST_UNIFORM(TestUniformDistributionNeg12Pos1DoubleSplit)

TEST_UNIFORM(TestUniformDistributionNeg12Pos1FullSplit)

TEST_NORMAL(TestNormalDistribution)

TEST_NORMAL(TestNormalDistributionWithSplit)

TEST_NORMAL(TestNormalDistributionWithDoubleSplit)

TEST_NORMAL(TestNormalDistributionWithFullSplit)

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithRandomInitialized()
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}