#include "GelfandProblem.hpp"

GelfandProblem *CreateGelfandProblem(const string &problemName) {
  if (problemName == "Gelfand") return new GelfandProblem();
  else Exit(problemName + " not found for SpaceDim=" + std::to_string(SpaceDimension))
}

std::shared_ptr<GelfandProblem> CreateGelfandProblemShared(const string &problemName) {
  return std::shared_ptr<GelfandProblem>(CreateGelfandProblem(problemName));
}
