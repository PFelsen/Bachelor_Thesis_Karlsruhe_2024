#ifndef BASICNONLINEARSOLVER_HPP
#define BASICNONLINEARSOLVER_HPP

#include "BasicSolver.hpp"
#include "Config.hpp"

class BasicFunction {
  double h;
public:
  BasicFunction(double numh = 1e-8) : h(numh) {}

  /// evaluates the function at point x
  virtual RVector operator()(const RVector &x) const = 0;

  /// evaluates the jacobian at point x (by default the jacobian is computed numerically)
  virtual RMatrix Jacobian(const RVector &x) const { return this->computeJacobian(x); }
private:
  /// numerical computation of the jacobian
  RMatrix computeJacobian(const RVector &x) const;
};

class BasicNewton {
  std::string linearSolverName;
  int verbose = 0;
  int iter;
  int max_iter = 10;
  int min_iter = 0;
  int LS_iter = 3;
  double Eps = 1e-8;
  double Red = 1e-8;
  double MinRed = 1e-9;
  double d;
public:
  BasicNewton(const std::string &linearSolverName) : linearSolverName(linearSolverName) {
    Config::Get("BasicNewtonVerbose", verbose);
    Config::Get("BasicNewtonSteps", max_iter);
    Config::Get("BasicNewtonLineSearchSteps", LS_iter);
    Config::Get("BasicNewtonEpsilon", Eps);
    Config::Get("BasicNewtonReduction", Red);
    Config::Get("BasicNewtonMinimalReduction", MinRed);
    Config::Get("BasicNewtonMinimalStepNumber", min_iter);
  }

  std::string Name() const { return "BasicNewtonMethod"; }

  void operator()(RVector &x, BasicFunction &f) const;
};


#endif // BASICNONLINEARSOLVER_HPP
