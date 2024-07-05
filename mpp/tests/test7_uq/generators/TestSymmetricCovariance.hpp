#ifndef TESTSYMMETRICCOVARIANCE_HPP
#define TESTSYMMETRICCOVARIANCE_HPP

#include "FFT.hpp"
#include "SymmetricCovariance.hpp"
#include "TestEnvironment.hpp"

template<typename T>
class TestCovariance : public Test {
protected:
  SymmetricCovariance<T> *covariance;

  explicit TestCovariance(SymmetricCovariance<T> *covariance) : covariance(covariance) {
    mpp_ba::SetTolerance(1e-6);
  }

  void TearDown() override { delete covariance; }
};

class TestCovariance1D : public TestCovariance<RVector> {
public:
  TestCovariance1D() : TestCovariance(new CovarianceFunction1D(1.0, 0.3, 1.0)) {}
};

class TestIsotropicCovariance2D : public TestCovariance<RMatrix> {
public:
  TestIsotropicCovariance2D() : TestCovariance(new CovarianceFunction2D(1.0, {0.3, 0.3}, 1.0, 1)) {}
};

class TestIsotropicCovariance3D : public TestCovariance<RTensor> {
public:
  TestIsotropicCovariance3D() :
      TestCovariance(new CovarianceFunction3D(1.0, {0.3, 0.3, 0.3}, 1.0, 1)) {}
};

#endif // TESTSYMMETRICCOVARIANCE_HPP
