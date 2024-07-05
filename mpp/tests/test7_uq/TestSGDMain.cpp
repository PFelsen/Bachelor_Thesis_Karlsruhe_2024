#include "TestSGDMain.hpp"

#include "StochasticGradientDescent.hpp"
#include "TestEnvironment.hpp"
#include "TestSGDMain.hpp"

class TestMainProgram : public TestWithParam<ConfigMap> {
protected:
  std::unique_ptr<StochasticGradientDescent> main;

  double referenceValue = 0.0;
  double referenceValueStepsize = 0.0;

  explicit TestMainProgram() {
    Config::Initialize(GetParam());
    main = std::make_unique<StochasticGradientDescent>();
    Config::Get("ReferenceValue", referenceValue);
    Config::Get("ReferenceValueStepsize", referenceValueStepsize);
  }

  void TearDown() override { Config::Close(); }
};

INSTANTIATE_TEST_SUITE_P(TestMainProgram, TestMainProgram, Values(DefaultSGDConfigMap));

TEST_P(TestMainProgram, TestRunMain) { EXPECT_NEAR(main->Compute(), referenceValue, 1e-6); }

TEST_P(TestMainProgram, TestValuesStepsize) {
  main->SetStepsize();
  EXPECT_NEAR(main->Stepsize(), referenceValueStepsize, 1e-6);
}

TEST_P(TestMainProgram, TestValuesIterationCounter) {
  double referenceValueIteration = main->Iteration();
  main->UpdateIteration();
  EXPECT_NEAR(main->Iteration(), referenceValueIteration + 1, 1e-6);
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithoutDefaultConfig()
                     .WithScreenLogging()
                     .WithRandomInitialized()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}