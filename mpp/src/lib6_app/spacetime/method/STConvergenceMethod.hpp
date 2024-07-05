#ifndef SPACETIME_STCONVERGENCEMAIN_HPP
#define SPACETIME_STCONVERGENCEMAIN_HPP

#include "STMethod.hpp"

class STConvergenceMethod : public STMethod {

public:
  STConvergenceMethod() : STConvergenceMethod(STMainBuilder().ReadConfig()){};

  STConvergenceMethod(STMainBuilder builder);

  void run() override;

};


#endif //SPACETIME_STCONVERGENCEMAIN_HPP
