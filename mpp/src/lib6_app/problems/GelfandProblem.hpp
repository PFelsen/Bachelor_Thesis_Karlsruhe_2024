#ifndef TUTORIAL_GELFANDPROBLEM_HPP
#define TUTORIAL_GELFANDPROBLEM_HPP

#include "IProblem.hpp"

class GelfandProblem : public IProblem {
public:
  double lambda = 1;
  double initlocation_left = -2.8;
  double initlocation_right = 3;
  double initvalue_left = 1;
  double initvalue_right = 1;
  GelfandProblem () : IProblem() {
    Config::Get("lambda", lambda);
    Config::Get("initvalue_left", initvalue_left);
    Config::Get("initvalue_right", initvalue_right);
    Config::Get("initlocation_left", initlocation_left);
    Config::Get("initlocation_right", initlocation_right);
  }
  string Name() const {
    return "GefandProblem with lambda = " + std::to_string(lambda);
  }
  double InitialValue (const Point& x) const {
    if (x[0] > 0) {
      double y = abs(initlocation_right - x[0]);
      if (y > 2) return 0;
      if (abs(x[1]) > 2) return 0;
      return 0.25 * initvalue_right * (2-y) * (2-abs(x[1]));
    } 
    double y = abs(initlocation_left - x[0]);
    if (y > 2) return 0;
    if (abs(x[1]) > 2) return 0;
    return 0.25 * initvalue_left * (2-y) * (2-abs(x[1]));
  }
  double DirichletValue (const Point& x) const { return 0; }
  double F (double u) const { return lambda * exp(u); }
  double DF (double u) const { return lambda * exp(u); }
};

GelfandProblem *CreateGelfandProblem(const string &problemName);

std::shared_ptr<GelfandProblem> CreateGelfandProblemShared(const std::string &problemName);

#endif // TUTORIAL_GELFANDPROBLEM_HPP
