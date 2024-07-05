
#include "TestVectorValuedMain.hpp"

INSTANTIATE_TEST_SUITE_P(Test, TestVectorValuedMainRates, ValuesIn(Parameters()));

TEST_VECTOR_POLYNOMIALS1D(TestVectorPolynomialsLagrange)
TEST_VECTOR_POLYNOMIALS1D(TestVectorPolynomialsDG)
TEST_VECTOR_POLYNOMIALS1D(TestVectorPolynomialsEG)

#if SpaceDimension >= 2
TEST_VECTOR_POLYNOMIALS2D(TestVectorPolynomialsLagrange)
TEST_VECTOR_POLYNOMIALS2D(TestVectorPolynomialsDG)
TEST_VECTOR_POLYNOMIALS2D(TestVectorPolynomialsEG)
#endif
#if SpaceDimension >= 3
TEST_VECTOR_POLYNOMIALS3D(TestVectorPolynomialsLagrange)
TEST_VECTOR_POLYNOMIALS3D(TestVectorPolynomialsDG)
TEST_VECTOR_POLYNOMIALS3D(TestVectorPolynomialsEG)
#endif

INITIALIZE_TEST_P(TestVectorPolynomialsLagrange)
INITIALIZE_TEST_P(TestVectorPolynomialsDG)
INITIALIZE_TEST_P(TestVectorPolynomialsEG)

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithoutDefaultConfig()
                     .WithConfPath(string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(string(ProjectMppDir) + "/conf/geo/")
                     .WithScreenLogging()
                     .WithFileLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
