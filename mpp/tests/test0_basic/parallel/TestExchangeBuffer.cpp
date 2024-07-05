#include "TestExchangeBuffer.hpp"
#include "parallel/Parallel.hpp"


EXCHANGEBUFFER_TESTS(TestExchangeBuffer);

EXCHANGEBUFFER_TESTS(TestExchangeBufferWithSplit)

EXCHANGEBUFFER_TESTS(TestExchangeBufferWithDoubleSplit)

EXCHANGEBUFFER_TESTS(TestExchangeBufferWithFullSplit)

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .
                 //            WithParallelListeners().
                 WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}