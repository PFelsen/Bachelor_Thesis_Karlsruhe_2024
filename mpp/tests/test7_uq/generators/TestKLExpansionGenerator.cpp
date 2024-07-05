#include "TestKLExpansionGenerator.hpp"

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithRandomInitialized()
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}