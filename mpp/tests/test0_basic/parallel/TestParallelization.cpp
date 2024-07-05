#include "TestParallelization.hpp"
#include "parallel/Parallel.hpp"


PPM_TESTS(TestPPM)

PPM_TESTS(TestPPMWithSplit)

PPM_TESTS(TestPPMWithDoubleSplit)

PPM_TESTS(TestPPMWithFullSplit)

TEST_F(TestPPMWithFullSplit, PrintInfo) {
  PPM->Barrier(0);
  for (int l = 0; l < PPM->CommsSize(); l++)
    PPM->PrintInfo(l);
  PPM->Barrier(0);
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithParallelListeners().WithScreenLogging().WithPPM())
      .RUN_ALL_MPP_TESTS();
}