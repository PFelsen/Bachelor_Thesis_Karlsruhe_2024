//
// Created by lstengel on 21.07.22.
//

#include "IVectorValuedAssemble.hpp"
#include "LagrangeVectorValuedAssemble.hpp"
#include "DGVectorValuedAssemble.hpp"
#include "EGVectorValuedAssemble.hpp"
#include "PDESolver.hpp"

IVectorValuedAssemble *CreateVectorValuedAssemble(IVectorValuedProblem &problem, const PDESolverConfig &conf) {

  if (conf.modelName.find("VectorValuedLagrange") != std::string::npos) {
    return new LagrangeVectorValuedAssemble(problem, conf.degree);
  }

  if (conf.modelName.find("EGVectorValuedAssemble") != std::string::npos) {
    return new EGVectorValuedAssemble(problem, conf.degree);
  }

  if (conf.modelName.find("DGVectorValuedAssemble") != std::string::npos) {
    return new DGVectorValuedAssemble(problem, conf.degree);
  }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<IVectorValuedAssemble> CreateVectorValuedAssembleUnique(
    IVectorValuedProblem &problem,  const PDESolverConfig &conf) {
  return std::unique_ptr<IVectorValuedAssemble>(CreateVectorValuedAssemble(problem, conf));
}

std::shared_ptr<IVectorValuedAssemble> CreateVectorValuedAssembleShared(
    IVectorValuedProblem &problem,  const PDESolverConfig &conf) {
  return std::shared_ptr<IVectorValuedAssemble>(CreateVectorValuedAssemble(problem, conf));
}