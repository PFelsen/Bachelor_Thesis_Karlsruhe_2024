#include "DGTransportAssemble.hpp"
#include "PDESolver.hpp"

ITransportAssemble *CreateTransportAssemble(
    const ITransportProblem &problem, const PDESolverConfig &conf) {

  if (conf.modelName.find("DG") != std::string::npos) {
    return new DGTransportAssemble(problem, conf.degree);
  }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<ITransportAssemble>
CreateTransportAssembleUnique(const ITransportProblem &problem, const PDESolverConfig &conf) {
  return std::unique_ptr<ITransportAssemble>(
      CreateTransportAssemble(problem, conf));
}

std::shared_ptr<ITransportAssemble> CreateTransportAssembleShared(
    const ITransportProblem &problem, const PDESolverConfig &conf) {
  return std::shared_ptr<ITransportAssemble>(
      CreateTransportAssemble(problem, conf));
}
