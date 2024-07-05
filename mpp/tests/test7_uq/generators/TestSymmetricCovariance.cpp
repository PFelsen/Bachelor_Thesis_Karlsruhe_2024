#include "TestSymmetricCovariance.hpp"

TEST_F(TestCovariance1D, TestToeplitzMatrix) {
  RVector toepRow(4);
  RVector toepCol(4);
  covariance->ToeplitzMatrix(toepRow, toepCol);
  RVector expectedToepRow{1.0, 0.3291930, 0.10836802, 0.0356740};
  RVector expectedToepCol{1.0, 0.3291930, 0.10836802, 0.0356740};
  EXPECT_EQ(expectedToepRow, toepRow);
  EXPECT_EQ(expectedToepCol, toepCol);

  // Has to be since covariance function is symmetric
  EXPECT_EQ(toepRow, toepCol);
}

TEST_F(TestCovariance1D, TestEmbeddedToeplitzMatrix) {
  RVector toepRow(4), toepCol(4);
  RVector circRow = covariance->EmbeddedToeplitzMatrix(toepRow, toepCol);

  RVector expectedCircRow{1.0, 0.3291930, 0.10836802, 0.0356740, 0.0356740, 0.10836802, 0.3291930};

  EXPECT_EQ(expectedCircRow, circRow);
}

TEST_F(TestCovariance1D, TestEigenValues) {
  RVector circRow{1.0, 0.1888756, 0.03567399, 0.03567399, 0.1888756};
  RVector eigenVal(circRow);
  FFT::RealToRealVector(circRow, eigenVal);
  RVector expectedEigenVal{1.44909919, 1.05900981, 0.7164406, 0.7164406, 1.05900981};
  EXPECT_EQ(expectedEigenVal, eigenVal);
}

TEST_F(TestIsotropicCovariance2D, TestToeplitzMatrix) {
  RMatrix toepRow(3, 3), toepCol(3, 3);
  covariance->ToeplitzMatrix(toepRow, toepCol);
  RMatrix expectedToepRow{{1.0, 0.1888756, 0.03567399},
                          {0.1888756, 0.03567399, 0.00673795},
                          {0.03567399, 0.00673795, 0.00127263}};
  RMatrix expectedToepCol{{1.0, 0.1888756, 0.03567399},
                          {0.1888756, 0.03567399, 0.00673795},
                          {0.03567399, 0.00673795, 0.00127263}};

  EXPECT_EQ(expectedToepRow, toepRow);
  EXPECT_EQ(expectedToepCol, toepCol);

  // Has to be since covariance function is symmetric
  EXPECT_EQ(toepRow, toepCol);

  // Has to be since isotropic (transpose)
  EXPECT_EQ(toepRow, toepRow);
  EXPECT_EQ(toepCol, toepCol);
}

class TestAnisotropicCovariance2D : public TestCovariance<RMatrix> {
public:
  TestAnisotropicCovariance2D() :
      TestCovariance(new CovarianceFunction2D(1.0, {0.3, 0.1}, 1.0, 1)) {}
};

TEST_F(TestAnisotropicCovariance2D, TestToeplitzMatrix) {
  RMatrix toepRow(3, 3), toepCol(3, 3);
  covariance->ToeplitzMatrix(toepRow, toepCol);

  // Has to be since covariance function is symmetric
  EXPECT_EQ(toepRow, toepCol);

  // Check for Anisotropy
  RMatrix testMatrix(toepRow);
  testMatrix -= toepRow.transpose();
  RVector testVector(testMatrix.diag());
  EXPECT_EQ(testVector, RVector(0.0, testVector.size()));
  //    RMatrix testMatrix2(-testMatrix);
  //    EXPECT_EQ(testMatrix, testMatrix2);
}

TEST_F(TestAnisotropicCovariance2D, TestEmbeddedToeplitzMatrix) {
  RMatrix toepRow(3, 3), toepCol(3, 3);
  RMatrix circRow = covariance->EmbeddedToeplitzMatrix(toepRow, toepCol);

  RMatrix expectedCircRow{{1.00000000e+00, 1.88875603e-01, 3.56739933e-02, 3.56739933e-02,
                           1.88875603e-01},
                          {6.73794700e-03, 1.27263380e-03, 2.40369476e-04, 2.40369476e-04,
                           1.27263380e-03},
                          {4.53999298e-05, 8.57493910e-06, 1.61959679e-06, 1.61959679e-06,
                           8.57493910e-06},
                          {4.53999298e-05, 8.57493910e-06, 1.61959679e-06, 1.61959679e-06,
                           8.57493910e-06},
                          {6.73794700e-03, 1.27263380e-03, 2.40369476e-04, 2.40369476e-04,
                           1.27263380e-03}};

  EXPECT_EQ(expectedCircRow, circRow);
}

TEST_F(TestAnisotropicCovariance2D, TestEigenValues) {
  RMatrix toepRow(3, 3), toepCol(3, 3);
  RMatrix circRow = covariance->EmbeddedToeplitzMatrix(toepRow, toepCol);

  RMatrix eigenVal(circRow);
  RMatrix expectedEV{{1.46875868, 1.07337707, 0.72616033, 0.72616033, 1.07337707},
                     {1.4550272, 1.06334203, 0.71937143, 0.71937143, 1.06334203},
                     {1.43334144, 1.04749396, 0.7086499, 0.7086499, 1.04749396},
                     {1.43334144, 1.04749396, 0.7086499, 0.7086499, 1.04749396},
                     {1.4550272, 1.06334203, 0.71937143, 0.71937143, 1.06334203}};
  FFT::RealToRealMatrix(circRow, eigenVal);
  EXPECT_EQ(eigenVal, expectedEV);
}

TEST_F(TestIsotropicCovariance3D, TestToeplitzMatrix) {
  RTensor toepRow(3, 3, 3), toepCol(3, 3, 3);
  covariance->ToeplitzMatrix(toepRow, toepCol);

  // Has to be since covariance function is symmetric
  EXPECT_EQ(toepRow, toepCol);
}

TEST_F(TestIsotropicCovariance3D, TestEmbeddedToeplitzMatrix) {
  RTensor toepRow(2), toepCol(2);
  for (int i = 0; i < toepRow.FirstDimension(); ++i) {
    for (int j = 0; j < toepRow.SecondDimension(); ++j) {
      for (int k = 0; k < toepRow.ThirdDimension(); ++k) {
        toepRow(i, j, k) = 4 - i - j - k;
        toepCol(i, j, k) = 4 - i - j - k;
      }
    }
  }
  RTensor circRow = BBMOCE(toepRow, toepCol);

  RTensor expectedCircRow{{{4, 3, 3}, {3, 2, 2}, {3, 2, 2}},
                          {{3, 2, 2}, {2, 1, 1}, {2, 1, 1}},
                          {{3, 2, 2}, {2, 1, 1}, {2, 1, 1}}};

  EXPECT_EQ(expectedCircRow, circRow);
}

TEST_F(TestIsotropicCovariance3D, TestEigenValues) {
  RTensor toepRow(2), toepCol(2);
  RTensor circRow = covariance->EmbeddedToeplitzMatrix(toepRow, toepCol);

  RTensor eigenVal(circRow);
  RTensor expectedEV{{{1.2296787651376919, 1.1068403803227072, 1.1068403803227072},
                      {1.1068403803227072, 0.99627289845550548, 0.99627289845550548},
                      {1.1068403803227072, 0.99627289845550548, 0.99627289845550548}},
                     {{1.1068403803227072, 0.99627289845550548, 0.99627289845550548},
                      {0.99627289845550548, 0.89675052143249967, 0.89675052143249967},
                      {0.99627289845550548, 0.89675052143249967, 0.89675052143249967}},
                     {{1.1068403803227072, 0.99627289845550548, 0.99627289845550548},
                      {0.99627289845550548, 0.89675052143249967, 0.89675052143249967},
                      {0.99627289845550548, 0.89675052143249967, 0.89675052143249967}}};
  FFT::RealToRealTensor(circRow, eigenVal);
  // Checked by python script. values are correct

  EXPECT_EQ(eigenVal, expectedEV);
}

/*
TEST_F(TestAnisotropicCovariance2D, TestEigenValues) {
    RMatrix toepRow(3, 3), toepCol(3, 3);
    RMatrix circRow = covariance->EmbeddedToeplitzMatrix(toepRow, toepCol);

    RMatrix eigenVal(circRow);
    RMatrix expectedEV{{1.46875868, 1.07337707, 0.72616033, 0.72616033, 1.07337707},
                       {1.4550272, 1.06334203, 0.71937143, 0.71937143, 1.06334203},
                       {1.43334144, 1.04749396, 0.7086499, 0.7086499, 1.04749396},
                       {1.43334144, 1.04749396, 0.7086499, 0.7086499, 1.04749396},
                       {1.4550272, 1.06334203, 0.71937143, 0.71937143, 1.06334203}};
    FFT::RealToRealMatrix(circRow, eigenVal);
    EXPECT_EQ(eigenVal, expectedEV);
}
*/

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithScreenLogging().WithPPM()).RUN_ALL_MPP_TESTS();
}