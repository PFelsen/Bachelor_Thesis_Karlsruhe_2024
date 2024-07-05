#include "PlateProblem.hpp"

PlateProblem *CreatePlateProblem(const string &problemName) {
  if (problemName == "Plate") return new PlateProblem();
  else Exit(problemName + " not found for SpaceDim=" + std::to_string(SpaceDimension))
}

std::shared_ptr<PlateProblem> CreatePlateProblemShared(const string &problemName) {
  return std::shared_ptr<PlateProblem>(CreatePlateProblem(problemName));
}
