#include "TestSparseGridGenerator.hpp"


SPARSEGRIDGENERATOR(TestClenshawCurtis1D)

SPARSEGRIDGENERATOR(TestClenshawCurtis2D)

SPARSEGRIDGENERATOR(TestClenshawCurtis5D)

SPARSEGRIDGENERATOR(TestClenshawCurtis1221D)

SPARSEGRIDGENERATOR(TestClenshawCurtis0212D)

SPARSEGRIDGENERATOR(TestLocalPolynomial1D)

SPARSEGRIDGENERATOR(TestLocalPolynomial2D)

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithScreenLogging().WithPPM()).RUN_ALL_MPP_TESTS();
}
