#include "AcousticProblems.hpp"
#include "DGAcousticAssemble.hpp"
#include "PDESolver.hpp"
IAcousticAssemble *
CreateAcousticAssemble(const AcousticProblem &problem, const PDESolverConfig &conf) {

  if (conf.modelName.find("DG") != std::string::npos) {
    return new DGAcousticAssemble(problem, conf.degree);
  }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<IAcousticAssemble>
CreateAcousticAssembleUnique(const AcousticProblem &problem, const PDESolverConfig &conf) {
  return std::unique_ptr<IAcousticAssemble>(
      CreateAcousticAssemble(problem, conf));
}

std::shared_ptr<IAcousticAssemble>
CreateAcousticAssembleShared(const AcousticProblem &problem, const PDESolverConfig &conf) {
  return std::shared_ptr<IAcousticAssemble>(
      CreateAcousticAssemble(problem, conf));
}