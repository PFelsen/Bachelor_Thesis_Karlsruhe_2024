#include "ErrorEstimator.hpp"


std::unique_ptr<ErrorEstimator>
CreateErrorEstimator(const std::string &name, STAssemble &assemble, LinearSolver &solver) {
  if (name == "ConformingResidual")
    return std::make_unique<L2ReliableResidualErrorEstimator>(assemble, solver);
  if (name == "DG"){
    return std::make_unique<DGErrorEstimator>(assemble, solver);
  }
  if (name == "Residual") return std::make_unique<ResidualErrorEstimator>(assemble, solver);
  if (name == "Dual") return std::make_unique<DualErrorEstimator>(assemble, solver);
  THROW("No Error Estimator know with name: " + name);
}