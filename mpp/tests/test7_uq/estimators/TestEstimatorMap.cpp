#include "EstimatorMap.hpp"
#include "TestEnvironment.hpp"

class TestEstimatorMap : public Test {
protected:
  std::map<int, WelfordAggregate> aggregateMap;

  EstimatorMap estimatorMap;

  TestEstimatorMap() :
      aggregateMap(
          {{3, WelfordAggregate(SampleCounter{1600, 0, 0},
                                CostInSeconds(WelfordData<double>(400, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                  1600, 0.0)),
                                QuantityOfInterest(WelfordData<double>(2.0e-02, 64.0e-02, 0.0, 0.0,
                                                                       0.0, 0.0, 0.0, 1600, 0.0)),
                                DeltaQuantityOfInterest(WelfordData<double>(2.0e-02, 64.0e-02, 0.0,
                                                                            0.0, 0.0, 0.0, 0.0,
                                                                            1600, 0.0)),
                                640000, 400)},
           {4, WelfordAggregate(SampleCounter{400, 0, 0},
                                CostInSeconds(WelfordData<double>(800, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                  400, 0.0)),
                                QuantityOfInterest(WelfordData<double>(1.0e-02, 16.0e-02, 0.0, 0.0,
                                                                       0.0, 0.0, 0.0, 400, 0.0)),
                                DeltaQuantityOfInterest(WelfordData<double>(1.0e-02, 16.0e-02, 0.0,
                                                                            0.0, 0.0, 0.0, 0.0, 400,
                                                                            0.0)),
                                320000, 800)},
           {5,
            WelfordAggregate(SampleCounter(100, 0, 0),
                             CostInSeconds(
                                 WelfordData<double>(1600, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100, 0.0)),
                             QuantityOfInterest(WelfordData<double>(0.5e-02, 4.0e-02, 0.0, 0.0, 0.0,
                                                                    0.0, 0.0, 100, 0.0)),
                             DeltaQuantityOfInterest(WelfordData<double>(0.5e-02, 4.0e-02, 0.0, 0.0,
                                                                         0.0, 0.0, 0.0, 100, 0.0)),
                             160000, 1600)},
           {6, WelfordAggregate(SampleCounter{25, 0, 0},
                                CostInSeconds(WelfordData<double>(3200, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                  0.0, 25, 0.0)),
                                QuantityOfInterest(WelfordData<double>(0.25e-02, 1.0e-02, 0.0, 0.0,
                                                                       0.0, 0.0, 0.0, 25, 0.0)),
                                DeltaQuantityOfInterest(WelfordData<double>(0.25e-02, 1.0e-02, 0.0,
                                                                            0.0, 0.0, 0.0, 0.0, 25,
                                                                            0.0)),
                                80000, 3200)}}),
      estimatorMap(EstimatorMap(aggregateMap, MLEstimatorConfig().WithTimeBudget(10000000))){};

  void TearDown() override { PPM->Barrier(0); }
};

TEST_F(TestEstimatorMap, TestExponentsAndErrors) {
  estimatorMap.UpdateExponentsAndErrors();
  EXPECT_FLOAT_EQ(estimatorMap.GetExponents().alpha, 1.0);
  EXPECT_FLOAT_EQ(estimatorMap.GetExponents().beta, 2.0);
  EXPECT_FLOAT_EQ(estimatorMap.GetExponents().gamma, 1.0);
  EXPECT_FLOAT_EQ(estimatorMap.GetErrors().disc, 0.25e-02);
  EXPECT_FLOAT_EQ(estimatorMap.GetErrors().input, 1.6e-3);
}

TEST_F(TestEstimatorMap, TestMonteCarloMapUpdateSampleCounter) {
  estimatorMap.UpdateSampleCounterOnMap(0.01, 0.5);
  for (auto &[level, estimator] : estimatorMap) {
    EXPECT_TRUE(estimator.SamplesToCompute() != 0);
    EXPECT_EQ(estimator.GetAggregate().ctr.M, aggregateMap.find(level)->second.ctr.M);
  }
}

TEST_F(TestEstimatorMap, TestMonteCarloMapAppendLevel) {
  int newLevel = 7;
  estimatorMap.UpdateExponentsAndErrors();
  estimatorMap.UpdateNewLevel(0.1, 0.5, newLevel);
  EXPECT_FLOAT_EQ(estimatorMap.find(newLevel)->second.SamplesToCompute(), 6); // floor(25 / 4)
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithScreenLogging().WithPPM()).RUN_ALL_MPP_TESTS();
}