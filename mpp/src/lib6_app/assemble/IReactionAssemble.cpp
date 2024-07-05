#include "DGReactionAssemble.hpp"
#include "PGReactionAssemble.hpp"
#include "PDESolver.hpp"

IReactionAssemble *CreateReactionAssemble(
    const IReactionProblem &problem, const PDESolverConfig &conf) {
  if (conf.modelName.find("PG") != std::string::npos) {
    return new PGReactionAssemble(problem, conf.degree);
  }

  if (conf.modelName.find("DG") != std::string::npos) {
    return new DGReactionAssemble(problem, conf.degree);
  }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<IReactionAssemble> CreateReactionAssembleUnique(
    const IReactionProblem &problem, const PDESolverConfig &conf) {
  return std::unique_ptr<IReactionAssemble>(
    CreateReactionAssemble(problem, conf));
}

std::shared_ptr<IReactionAssemble> CreateReactionAssembleShared(
    const IReactionProblem &problem, const PDESolverConfig &conf) {
  return std::shared_ptr<IReactionAssemble>(
    CreateReactionAssemble(problem, conf));
}