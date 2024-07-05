#ifndef SPACETIME_STADAPTIVEMAIN_HPP
#define SPACETIME_STADAPTIVEMAIN_HPP


#include "STMethod.hpp"
#include "ErrorEstimator.hpp"

class STAdaptiveMethod : public STMethod {
protected:

  int refinement;

  std::string problemName;
  bool plot;

public:

  STAdaptiveMethod();

  STAdaptiveMethod(STMainBuilder builder);

  void run() override;

  void handleAdaptivity(const Vector &Eta) {
    mout.StartBlock("Adaptivity");
    double theta = 0;
    double theta_min = 0;
    double theta_factor = 1.0;
    std::string refine_by;

    Config::Get("theta", theta);
    Config::Get("theta_factor", theta_factor);
    Config::Get("theta_min", theta_min);
    Config::Get("refine_by", refine_by);

    theta *= pow(theta_factor, Eta.Level().adaptivityLevel);


    assemble->AdaptPolynomialDegrees(Eta, theta, theta_min, refine_by);
    assemble->AdaptQuadrature(Eta.Level().NextInAdaptivity(), -1, true);
    mout.EndBlock();
  }


};


#endif //SPACETIME_STADAPTIVEMAIN_HPP
