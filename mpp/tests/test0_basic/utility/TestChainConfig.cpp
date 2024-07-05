#include "TestEnvironment.hpp"


#include "gtest/gtest.h"
#include "utility/Config.hpp"

using namespace std;
using namespace ::testing;

struct ChainParams {
public:
  std::string value;
  std::string resultingValue;
};

class ChainConfigTest : public TestWithParam<ChainParams> {
public:
  ChainConfigTest() {}

  void SetUp() override {}

  void TearDown() override { Config::Close(); }
};

TEST_P(ChainConfigTest, ConfigTest) {
  std::string arg = "A=" + GetParam().value;
  int *argc = new int{2};
  char *argv[] = {(char *)"TestExecutable", (char *)arg.c_str(), nullptr};

  Config::AddChainConfig({"A", "B"}, {"key", "1"});
  Config::AddChainConfig({"A", "C"}, {"key", "2"});
  Config::AddChainConfig({"A", "D"}, {"key", "3"});

  Config::Initialize(argc, argv);

  string value = "N//A";
  Config::Get("key", value);
  EXPECT_EQ(value, GetParam().resultingValue);
}

INSTANTIATE_TEST_SUITE_P(ConfigTest, ChainConfigTest,
                         Values(ChainParams{"B", "1"}, ChainParams{"A", "N//A"}));

int main(int argc, char **argv) {
  Config::SetConfPath(Config::GetSourcePath() + "/conf/");

  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
