#include "TestRestriction.hpp"

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging().WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
