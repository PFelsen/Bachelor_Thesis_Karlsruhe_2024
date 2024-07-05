#include "TestEnvironment.hpp"

#include "TestSpectrum.hpp"

TEST_P(RealSpectrum_A_Test, EVRealTest) {
  Eigenvalues lambda;

  EVreal(A, lambda);
  EXPECT_EQ(lambda, expectedEigenvalues);

  //    EVreal(A_r, lambda);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
}

TEST_P(RealSpectrum_A_Test, EVRealBothTest) {
  Eigenvalues lambda;
  Eigenvectors e;

  EVreal(A, lambda, e);
  EXPECT_EQ(lambda, expectedEigenvalues);
  EXPECT_EQ(e, expectedEigenvectors);

  //    EVreal(A_r, lambda, e);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
  //    EXPECT_EQ(e, expectedEigenvectors);
}

TEST_P(RealSpectrum_A_Test, SpectrumTest) {
  Spectrum spec(A);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

TEST_P(RealSpectrum_A_Test, SpectrumBothTest) {
  Spectrum spec(A, true);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_EQ(spec.getEigenvectors(), expectedEigenvectors);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

INSTANTIATE_TEST_SUITE_P(SpectrumTest, RealSpectrum_A_Test, ValuesIn(RealSpectrum_A_TestCases));

TEST_P(RealSpectrum_AB_Test, EVRealTest) {
  Eigenvalues lambda;

  EVreal(A, B, lambda);
  EXPECT_EQ(lambda, expectedEigenvalues);

  //    EVreal(A_r, B_r, lambda);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
}

TEST_P(RealSpectrum_AB_Test, EVRealBothTest) {
  Eigenvalues lambda;
  Eigenvectors e;

  EVreal(A, B, lambda, e);
  EXPECT_EQ(lambda, expectedEigenvalues);
  EXPECT_EQ(e, expectedEigenvectors);

  //    EVreal(A_r, B_r, lambda, e);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
  //    EXPECT_EQ(e, expectedEigenvectors);
}

TEST_P(RealSpectrum_AB_Test, SpectrumTest) {
  Spectrum spec(A, B);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

TEST_P(RealSpectrum_AB_Test, SpectrumBothTest) {
  Spectrum spec(A, B, true);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_EQ(spec.getEigenvectors(), expectedEigenvectors);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

INSTANTIATE_TEST_SUITE_P(SpectrumTest, RealSpectrum_AB_Test, ValuesIn(RealSpectrum_AB_TestCases));

TEST_P(ComplexSpectrum_A_Test, EVComplexTest) {
  Eigenvalues lambda;

  EVcomplex(A, lambda);
  EXPECT_EQ(lambda, expectedEigenvalues);

  //    EVcomplex(A_c, lambda);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
}

TEST_P(ComplexSpectrum_A_Test, EVComplexBothTest) {
  Eigenvalues lambda;
  CEigenvectors e;

  EVcomplex(A, lambda, e);
  EXPECT_EQ(lambda, expectedEigenvalues);
  EXPECT_EQ(e, expectedEigenvectors);

  //    EVcomplex(A_c, lambda, e);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
  //    EXPECT_EQ(e, expectedEigenvectors);
}

TEST_P(ComplexSpectrum_A_Test, SpectrumTest) {
  Spectrum spec(A);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

TEST_P(ComplexSpectrum_A_Test, SpectrumBothTest) {
  Spectrum spec(A, true);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_EQ(spec.getCEigenvectors(), expectedEigenvectors);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

INSTANTIATE_TEST_SUITE_P(SpectrumTest, ComplexSpectrum_A_Test,
                         ValuesIn(ComplexSpectrum_A_TestCases));

TEST_P(ComplexSpectrum_AB_Test, EVComplexTest) {
  Eigenvalues lambda;

  EVcomplex(A, B, lambda);
  EXPECT_EQ(lambda, expectedEigenvalues);

  //    EVcomplex(A_c, B_c, lambda);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
}

TEST_P(ComplexSpectrum_AB_Test, EVComplexBothTest) {
  Eigenvalues lambda;
  CEigenvectors e;

  EVcomplex(A, B, lambda, e);
  EXPECT_EQ(lambda, expectedEigenvalues);
  EXPECT_EQ(e, expectedEigenvectors);

  //    EVcomplex(A_c, B_c, lambda, e);
  //    EXPECT_EQ(lambda, expectedEigenvalues);
  //    EXPECT_EQ(e, expectedEigenvectors);
}

TEST_P(ComplexSpectrum_AB_Test, SpectrumTest) {
  Spectrum spec(A, B);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

TEST_P(ComplexSpectrum_AB_Test, SpectrumBothTest) {
  Spectrum spec(A, B, true);
  EXPECT_EQ(spec(), expectedEigenvalues);
  EXPECT_EQ(spec.getCEigenvectors(), expectedEigenvectors);
  EXPECT_NEAR(spec.min(), expectedMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.max(), expectedMaxEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmin(), expectedAbsMinEigenvalue, valueTol);
  EXPECT_NEAR(spec.absmax(), expectedAbsMaxEigenvalue, valueTol);
}

INSTANTIATE_TEST_SUITE_P(SpectrumTest, ComplexSpectrum_AB_Test,
                         ValuesIn(ComplexSpectrum_AB_TestCases));

int main(int argc, char **argv) {
  combineTestCases();
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
