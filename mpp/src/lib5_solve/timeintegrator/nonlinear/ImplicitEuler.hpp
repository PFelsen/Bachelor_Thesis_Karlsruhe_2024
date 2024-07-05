#ifndef _EULER_H_
#define _EULER_H_

#include "Newton.hpp"
#include "TimeIntegrator.hpp"
#include "TimeSeries.hpp"

class ImplicitEuler : public NonLinearTimeIntegrator {
  int extrapolate = 1;

  int max_estep = 100000;
public:
  ImplicitEuler(NonLinearSolver *nonLinearSolver) : NonLinearTimeIntegrator(nonLinearSolver) {
    Config::Get("EulerVerbose", verbose);
    Config::Get("EulerExtrapolate", extrapolate);
    Config::Get("EulerSteps", max_estep);
  }

  int GetOrder() const override { return 1; };

  const char *Name() const override { return "ImplicitEuler"; };

  bool Method(INonLinearTimeAssemble &assemble, Vector &u) override;
};

#endif // _EULER_H_
