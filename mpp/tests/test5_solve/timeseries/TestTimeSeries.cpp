#include "TestTimeSeries.hpp"
#include "TestEnvironment.hpp"
#include "TimeSeries.hpp"

class TestTimeSeries : public Test {
protected:
  std::unique_ptr<TimeSeries> ts;

  void SetUp() override {
    ConfigMap timeSeriesConfig = testConfigMap;
    Config::Initialize(timeSeriesConfig);
    ts = std::make_unique<TimeSeries>();
  }

  void TearDown() override { Config::Close(); }
};

class TestTimeSeries2 : public Test {
protected:
  std::unique_ptr<TimeSeries> ts;

  void SetUp() override { ts = std::make_unique<TimeSeries>(0, 0.1, 0.1 * std::pow(2, -3)); }

  void TearDown() override { Config::Close(); }

  int GetTrueSteps() {
    int count = 0;
    while (!(*ts).IsFinished()) {
      ts->NextTimeStep();
      count++;
    }
    return count;
  }

  double GetEndTime() {
    while (!(*ts).IsFinished()) {
      ts->NextTimeStep();
    }
    return ts->Time();
  }
};

TEST_F(TestTimeSeries, TestConfigConstructor) {
  ASSERT_EQ(ts->FirstTStep(), 1.0);
  ASSERT_EQ(ts->LastTStep(), 3.0);
  ASSERT_EQ(ts->StepSize(), 2.0);
}

TEST_F(TestTimeSeries2, TestEndTime) { ASSERT_EQ(ts->LastTStep(), GetEndTime()); }

TEST_F(TestTimeSeries2, TestSteps) { ASSERT_EQ(ts->Steps(), GetTrueSteps()); }

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
