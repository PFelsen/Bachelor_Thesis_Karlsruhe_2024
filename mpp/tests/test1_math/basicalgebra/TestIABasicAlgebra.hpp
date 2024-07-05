#ifndef TESTIABASICALGEBRA_HPP
#define TESTIABASICALGEBRA_HPP

#include "RVector.hpp"
#include "SaveLoad.hpp"
#include "TestEnvironment.hpp"

#include <ctime>
#include <utility>
#include <vector>

const double tol = 1e-14;

using std::vector;

struct IABasicAlgebraVectorTestParameter {
  int dim;
};

class IABasicAlgebraVectorTest : public TestWithParam<IABasicAlgebraVectorTestParameter> {
protected:
  int dim;
public:
  IABasicAlgebraVectorTest() {
    dim = GetParam().dim;
    mpp_ba::SetTolerance(1e-14);
  }
};

struct IABasicAlgebraMatrixTestParameter {
  int rows;
  int cols;
};

class IABasicAlgebraMatrixTest : public TestWithParam<IABasicAlgebraMatrixTestParameter> {
protected:
  int rows;
  int cols;
public:
  IABasicAlgebraMatrixTest() {
    rows = GetParam().rows;
    cols = GetParam().cols;
    mpp_ba::SetTolerance(1e-14);
  }
};

#endif // TESTIABASICALGEBRA_HPP
