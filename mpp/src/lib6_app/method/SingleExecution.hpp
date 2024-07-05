#ifndef MPP_SINGLEEXECUTION_HPP
#define MPP_SINGLEEXECUTION_HPP

#include "PDESolver.hpp"

class SingleExecution {
protected:
  PDESolverConfig conf;

  std::map<std::string, double> results;

public:
  explicit SingleExecution(PDESolverConfig conf) : conf(std::move(conf)) {}

  std::map<std::string, double> GetResults() const { return results; }

  virtual void Method();
};

#endif //MPP_SINGLEEXECUTION_HPP
