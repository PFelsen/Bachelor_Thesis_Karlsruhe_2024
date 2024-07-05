#ifndef TESTSPECTRUM_H
#define TESTSPECTRUM_H

#include "gtest/gtest.h"

#include "spectrum/Spectrum.hpp"

#include <utility>
#include <math.h>


using std::vector;
using namespace ::testing;

const double valueTol = 1e-10;

struct RealSpectrum_A_TestParameter {
  int dim;
  vector<vector<double>> A;
  vector<double> expectedEigenvalues;
  vector<vector<double>> expectedEigenvectors;
  double expectedMinEigenvalue;
  double expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue;
  double expectedAbsMaxEigenvalue;
};

static vector<RealSpectrum_A_TestParameter>
    RealSpectrum_A_TestCases{RealSpectrum_A_TestParameter{2,
                                                          {{3.0, 2.0}, {2.0, 6.0}},
                                                          {2.0, 7.0},
                                                          {{-2.0 / std::sqrt(5.0),
                                                            1.0 / std::sqrt(5.0)},
                                                           {1.0 / std::sqrt(5.0),
                                                            2.0 / std::sqrt(5.0)}},
                                                          2.0,
                                                          7.0,
                                                          2.0,
                                                          7.0},
                             RealSpectrum_A_TestParameter{2,
                                                          {{1.0, -1.0}, {-1.0, 1.0}},
                                                          {0.0, 2.0},
                                                          {{-1.0 / std::sqrt(2.0),
                                                            -1.0 / std::sqrt(2.0)},
                                                           {-1.0 / std::sqrt(2.0),
                                                            1.0 / std::sqrt(2.0)}},
                                                          0.0,
                                                          2.0,
                                                          0.0,
                                                          2.0}};

class RealSpectrum_A_Test : public TestWithParam<RealSpectrum_A_TestParameter> {
protected:
  int dim;
  SymRMatrix A;
  RMatrix A_r;

  Eigenvalues expectedEigenvalues;
  RMatrix expectedEigenvectors;
  double expectedMinEigenvalue, expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue, expectedAbsMaxEigenvalue;

  RealSpectrum_A_Test() {
    dim = GetParam().dim;
    A.resize(dim);
    A_r.resize(dim);
    expectedEigenvalues.resize(dim);
    expectedEigenvectors.resize(dim);
    for (int n = 0; n < dim; ++n) {
      expectedEigenvalues[n] = GetParam().expectedEigenvalues[n];
      for (int m = 0; m < dim; ++m) {
        A(n, m) = GetParam().A[n][m];
        A_r[n][m] = GetParam().A[n][m];
        expectedEigenvectors[n][m] = GetParam().expectedEigenvectors[n][m];
      }
    }
    expectedMinEigenvalue = GetParam().expectedMinEigenvalue;
    expectedMaxEigenvalue = GetParam().expectedMaxEigenvalue;
    expectedAbsMinEigenvalue = GetParam().expectedAbsMinEigenvalue;
    expectedAbsMaxEigenvalue = GetParam().expectedAbsMaxEigenvalue;
  }
};

struct RealSpectrum_AB_TestParameter {
  int dim;
  vector<vector<double>> A;
  vector<vector<double>> B;
  vector<double> expectedEigenvalues;
  vector<vector<double>> expectedEigenvectors;
  double expectedMinEigenvalue;
  double expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue;
  double expectedAbsMaxEigenvalue;
};

static vector<RealSpectrum_AB_TestParameter> RealSpectrum_AB_TestCases{
    // TODO: Add test cases
};

class RealSpectrum_AB_Test : public TestWithParam<RealSpectrum_AB_TestParameter> {
protected:
  int dim;
  SymRMatrix A, B;
  RMatrix A_r, B_r;

  Eigenvalues expectedEigenvalues;
  RMatrix expectedEigenvectors;
  double expectedMinEigenvalue, expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue, expectedAbsMaxEigenvalue;

  RealSpectrum_AB_Test() {
    dim = GetParam().dim;
    A.resize(dim);
    B.resize(dim);
    A_r.resize(dim);
    B_r.resize(dim);
    expectedEigenvalues.resize(dim);
    expectedEigenvectors.resize(dim);
    for (int n = 0; n < dim; ++n) {
      expectedEigenvalues[n] = GetParam().expectedEigenvalues[n];
      for (int m = 0; m < dim; ++m) {
        A(n, m) = GetParam().A[n][m];
        A_r[n][m] = GetParam().A[n][m];
        B(n, m) = GetParam().B[n][m];
        B_r[n][m] = GetParam().B[n][m];
        expectedEigenvectors[n][m] = GetParam().expectedEigenvectors[n][m];
      }
    }
    expectedMinEigenvalue = GetParam().expectedMinEigenvalue;
    expectedMaxEigenvalue = GetParam().expectedMaxEigenvalue;
    expectedAbsMinEigenvalue = GetParam().expectedAbsMinEigenvalue;
    expectedAbsMaxEigenvalue = GetParam().expectedAbsMaxEigenvalue;
  }
};

struct ComplexSpectrum_A_TestParameter {
  int dim;
  vector<vector<std::complex<double>>> A;
  vector<double> expectedEigenvalues;
  vector<vector<std::complex<double>>> expectedEigenvectors;
  double expectedMinEigenvalue;
  double expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue;
  double expectedAbsMaxEigenvalue;
};

static vector<ComplexSpectrum_A_TestParameter> ComplexSpectrum_A_TestCases{
    // TODO: Add test cases
};

class ComplexSpectrum_A_Test : public TestWithParam<ComplexSpectrum_A_TestParameter> {
protected:
  int dim;
  HermCMatrix A;
  CMatrix A_c;

  Eigenvalues expectedEigenvalues;
  CMatrix expectedEigenvectors;
  double expectedMinEigenvalue, expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue, expectedAbsMaxEigenvalue;

  ComplexSpectrum_A_Test() {
    dim = GetParam().dim;
    A.resize(dim);
    A_c.resize(dim);
    expectedEigenvalues.resize(dim);
    expectedEigenvectors.resize(dim);
    for (int n = 0; n < dim; ++n) {
      expectedEigenvalues[n] = GetParam().expectedEigenvalues[n];
      for (int m = 0; m < dim; ++m) {

        A(GetParam().A[n][m], n, m);
        A_c[n][m] = GetParam().A[n][m];
        expectedEigenvectors[n][m] = GetParam().expectedEigenvectors[n][m];
      }
    }
    expectedMinEigenvalue = GetParam().expectedMinEigenvalue;
    expectedMaxEigenvalue = GetParam().expectedMaxEigenvalue;
    expectedAbsMinEigenvalue = GetParam().expectedAbsMinEigenvalue;
    expectedAbsMaxEigenvalue = GetParam().expectedAbsMaxEigenvalue;
  }
};

struct ComplexSpectrum_AB_TestParameter {
  int dim;
  vector<vector<std::complex<double>>> A;
  vector<vector<std::complex<double>>> B;
  vector<double> expectedEigenvalues;
  vector<vector<std::complex<double>>> expectedEigenvectors;
  double expectedMinEigenvalue;
  double expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue;
  double expectedAbsMaxEigenvalue;
};

static vector<ComplexSpectrum_AB_TestParameter> ComplexSpectrum_AB_TestCases{
    // TODO: Add test cases
};

class ComplexSpectrum_AB_Test : public TestWithParam<ComplexSpectrum_AB_TestParameter> {
protected:
  int dim;
  HermCMatrix A, B;
  CMatrix A_c, B_c;

  Eigenvalues expectedEigenvalues;
  CMatrix expectedEigenvectors;
  double expectedMinEigenvalue, expectedMaxEigenvalue;
  double expectedAbsMinEigenvalue, expectedAbsMaxEigenvalue;

  ComplexSpectrum_AB_Test() {
    dim = GetParam().dim;
    A.resize(dim);
    A_c.resize(dim);
    B.resize(dim);
    B_c.resize(dim);
    expectedEigenvalues.resize(dim);
    expectedEigenvectors.resize(dim);
    for (int n = 0; n < dim; ++n) {
      expectedEigenvalues[n] = GetParam().expectedEigenvalues[n];
      for (int m = 0; m < dim; ++m) {
        A(GetParam().A[n][m], n, m);
        A_c[n][m] = GetParam().A[n][m];
        B(GetParam().B[n][m], n, m);
        B_c[n][m] = GetParam().B[n][m];
        expectedEigenvectors[n][m] = GetParam().expectedEigenvectors[n][m];
      }
    }
    expectedMinEigenvalue = GetParam().expectedMinEigenvalue;
    expectedMaxEigenvalue = GetParam().expectedMaxEigenvalue;
    expectedAbsMinEigenvalue = GetParam().expectedAbsMinEigenvalue;
    expectedAbsMaxEigenvalue = GetParam().expectedAbsMaxEigenvalue;
  }
};

void combineTestCases() {
  for (int i = 0; i < RealSpectrum_A_TestCases.size(); ++i) {
    RealSpectrum_AB_TestCases.push_back(
        RealSpectrum_AB_TestParameter{RealSpectrum_A_TestCases[i].dim,
                                      RealSpectrum_A_TestCases[i].A,
                                      {{1.0, 0.0}, {0.0, 1.0}},
                                      RealSpectrum_A_TestCases[i].expectedEigenvalues,
                                      RealSpectrum_A_TestCases[i].expectedEigenvectors,
                                      RealSpectrum_A_TestCases[i].expectedMinEigenvalue,
                                      RealSpectrum_A_TestCases[i].expectedMaxEigenvalue,
                                      RealSpectrum_A_TestCases[i].expectedAbsMinEigenvalue,
                                      RealSpectrum_A_TestCases[i].expectedAbsMaxEigenvalue});

    vector<vector<std::complex<double>>> A_complex(RealSpectrum_A_TestCases[i].A.size());
    for (int k = 0; k < A_complex.size(); ++k) {
      A_complex[k].resize(RealSpectrum_A_TestCases[i].A[k].size());
      for (int l = 0; l < A_complex[k].size(); ++l)
        A_complex[k][l] = std::complex<double>{RealSpectrum_A_TestCases[i].A[k][l], 0.0};
    }

    vector<vector<std::complex<double>>> e_complex(
        RealSpectrum_A_TestCases[i].expectedEigenvectors.size());
    for (int k = 0; k < e_complex.size(); ++k) {
      e_complex[k].resize(RealSpectrum_A_TestCases[i].expectedEigenvectors[k].size());
      for (int l = 0; l < e_complex[k].size(); ++l)
        e_complex[k][l] = RealSpectrum_A_TestCases[i].expectedEigenvectors[k][l];
    }

    ComplexSpectrum_A_TestCases.push_back(
        ComplexSpectrum_A_TestParameter{RealSpectrum_A_TestCases[i].dim, A_complex,
                                        RealSpectrum_A_TestCases[i].expectedEigenvalues, e_complex,
                                        RealSpectrum_A_TestCases[i].expectedMinEigenvalue,
                                        RealSpectrum_A_TestCases[i].expectedMaxEigenvalue,
                                        RealSpectrum_A_TestCases[i].expectedAbsMinEigenvalue,
                                        RealSpectrum_A_TestCases[i].expectedAbsMaxEigenvalue});
  }

  for (int i = 0; i < ComplexSpectrum_A_TestCases.size(); ++i) {
    ComplexSpectrum_AB_TestCases.push_back(
        ComplexSpectrum_AB_TestParameter{ComplexSpectrum_A_TestCases[i].dim,
                                         ComplexSpectrum_A_TestCases[i].A,
                                         {{{1.0, 0.0}, {0.0, 0.0}}, {{0.0, 0.0}, {1.0, 0.0}}},
                                         ComplexSpectrum_A_TestCases[i].expectedEigenvalues,
                                         ComplexSpectrum_A_TestCases[i].expectedEigenvectors,
                                         ComplexSpectrum_A_TestCases[i].expectedMinEigenvalue,
                                         ComplexSpectrum_A_TestCases[i].expectedMaxEigenvalue,
                                         ComplexSpectrum_A_TestCases[i].expectedAbsMinEigenvalue,
                                         ComplexSpectrum_A_TestCases[i].expectedAbsMaxEigenvalue});
  }
}

#endif // TESTSPECTRUM_H
