#include "LagrangeEllipticAssemble.hpp"
#include "HybridEllipticAssemble.hpp"
#include "MixedEllipticAssemble.hpp"
#include "DGEllipticAssemble.hpp"
#include "EGEllipticAssemble.hpp"
#include "PDESolver.hpp"


IEllipticAssemble *
CreateEllipticAssemble(const IEllipticProblem &problem, const PDESolverConfig &conf) {
  if (conf.modelName.find("Lagrange") != std::string::npos) {

    return new LagrangeEllipticAssemble(problem, conf.degree);
  }

  if (conf.modelName.find("Mixed") != std::string::npos) {
    return new MixedEllipticAssemble(problem);
  }

  if (conf.modelName.find("Hybrid") != std::string::npos) {
    return new HybridEllipticAssemble(problem);
  }

  if (conf.modelName.find("DG") != std::string::npos) {
    return new DGEllipticAssemble(problem, conf.degree);
  }

  if (conf.modelName.find("EG") != std::string::npos) {
    return new EGEllipticAssemble(problem, conf.degree);
  }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<IEllipticAssemble>
CreateEllipticAssembleUnique(const IEllipticProblem &problem, const PDESolverConfig &conf) {
  return std::unique_ptr<IEllipticAssemble>( CreateEllipticAssemble(problem, conf));
}

std::shared_ptr<IEllipticAssemble> CreateEllipticAssembleShared(
    const IEllipticProblem &problem, const PDESolverConfig &conf) {
  return std::shared_ptr<IEllipticAssemble>(CreateEllipticAssemble(problem, conf));
}