#ifndef SPACETIME_STADAPTIVECONVERGENCEMAIN_HPP
#define SPACETIME_STADAPTIVECONVERGENCEMAIN_HPP


#include "STMethod.hpp"
#include "ErrorEstimator.hpp"

class STAdaptiveConvergenceMethod : public STMethod {
protected:

  int refinement;

  std::string problemName;
  bool plot;

  int currentRun = 0;

public:

  STAdaptiveConvergenceMethod();

  STAdaptiveConvergenceMethod(STMainBuilder builder);

  void run() override;

  void handleAdaptivity(const Vector &Eta);


};


#endif //SPACETIME_STADAPTIVECONVERGENCEMAIN_HPP
