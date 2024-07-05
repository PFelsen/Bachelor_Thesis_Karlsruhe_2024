#ifndef TESTBASICALGEBRA_HPP
#define TESTBASICALGEBRA_HPP

#include "TestEnvironment.hpp"
#include "basicalgebra/RVector.hpp"
#include "utility/SaveLoad.hpp"

#include <complex>
#include <ctime>
#include <utility>
#include <vector>

const double tol = 1e-14;

using std::vector;

struct BasicAlgebraVectorTestParameter {
  int dim;
};

class BasicAlgebraVectorTest : public TestWithParam<BasicAlgebraVectorTestParameter> {
protected:
  int dim;
public:
  BasicAlgebraVectorTest() {
    dim = GetParam().dim;
    mpp_ba::SetTolerance(1e-14);
  }
};

struct BasicAlgebraMatrixTestParameter {
  int rows;
  int cols;
};

class BasicAlgebraMatrixTest : public TestWithParam<BasicAlgebraMatrixTestParameter> {
protected:
  int rows;
  int cols;
public:
  BasicAlgebraMatrixTest() {
    rows = GetParam().rows;
    cols = GetParam().cols;
    mpp_ba::SetTolerance(1e-14);
  }
};

struct BasicAlgebraTensorTestParameter {
  int FirstComponents;
  int SecondComponents;
  int ThirdComponents;
};

class BasicAlgebraTensorTest : public TestWithParam<BasicAlgebraTensorTestParameter> {
protected:
  int FirstComponents;
  int SecondComponents;
  int ThirdComponents;
public:
  BasicAlgebraTensorTest() {
    FirstComponents = GetParam().FirstComponents;
    SecondComponents = GetParam().SecondComponents;
    ThirdComponents = GetParam().ThirdComponents;
    mpp_ba::SetTolerance(1e-14);
  }
};

#endif // TESTBASICALGEBRA_HPP
