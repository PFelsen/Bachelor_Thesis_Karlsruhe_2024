#ifndef EIGENSOLVERCREATOR_HPP
#define EIGENSOLVERCREATOR_HPP

#include "LOBPCG.hpp"
#include "LOBPCGSelective.hpp"
#include "RitzGalerkin.hpp"

class EigenSolverCreator {
  std::string name;
  std::string linearSolver = "GMRES";
  std::string preconditioner = "PS";
  std::string linearSolver2 = "GMRES";
  std::string preconditioner2 = "PS";
  std::string defectType = "Default";
  int maxstep = 30;
  int printSteps = 1;
  int verbose = 3;
  double eps = 1e-5;
public:
  EigenSolverCreator(std::string name, std::string prefix = "") : name(name) {
    Config::Get(prefix + "ESolverLinearSolver", linearSolver);
    Config::Get(prefix + "ESolverPreconditioner", preconditioner);
    Config::Get(prefix + "ESolverLinearSolver2", linearSolver2);
    Config::Get(prefix + "ESolverPreconditioner2", preconditioner2);
    Config::Get(prefix + "ESolverMaxStep", maxstep);
    Config::Get(prefix + "ESolverPrintSteps", printSteps);
    Config::Get(prefix + "ESolverEpsilon", eps);
    Config::Get(prefix + "ESolverDefectType", defectType);
    Config::Get(prefix + "ESolverVerbose", verbose);
  }

  EigenSolverCreator &WithLinearSolver(std::string linearSolver) {
    this->linearSolver = linearSolver;
    return *this;
  }

  EigenSolverCreator &WithPreconditioner(std::string preconditioner) {
    this->preconditioner = preconditioner;
    return *this;
  }

  EigenSolverCreator &WithLinearSolver2(std::string linearSolver) {
    this->linearSolver2 = linearSolver;
    return *this;
  }

  EigenSolverCreator &WithPreconditioner2(std::string preconditioner) {
    this->preconditioner2 = preconditioner;
    return *this;
  }

  EigenSolverCreator &WithVerbose(int verbose) {
    this->verbose = verbose;
    return *this;
  }

  EigenSolverCreator &WithMaxStep(int maxstep) {
    this->maxstep = maxstep;
    return *this;
  }

  EigenSolverCreator &WithPrintSteps(int printSteps) {
    this->printSteps = printSteps;
    return *this;
  }

  EigenSolverCreator &WithEpsilon(double eps) {
    this->eps = eps;
    return *this;
  }

  EigenSolverCreator &WithDefectType(std::string type) {
    this->defectType = type;
    return *this;
  }

  IEigenSolver *Create() {
    std::unique_ptr<LinearSolver> solver(
        GetLinearSolver(linearSolver, GetPC(preconditioner), "ESolverLinear"));
    std::unique_ptr<LinearSolver> solver2 = nullptr;
    if (defectType == "OperatorB")
      solver2 = std::unique_ptr<LinearSolver>(
          GetLinearSolver(linearSolver2, GetPC(preconditioner2), "ESolverLinear"));

    if (name == "RitzGalerkin")
      return new RitzGalerkin(std::move(solver), std::move(solver2), verbose, eps, maxstep,
                              printSteps);
    if (name == "LOBPCG")
      return new LOBPCG(std::move(solver), std::move(solver2), verbose, eps, maxstep, printSteps);
    if (name == "LOBPCGExtended")
      return new LOBPCGExtended(std::move(solver), std::move(solver2), verbose, eps, maxstep,
                                printSteps);
    if (name == "LOBPCGSelective")
      return new LOBPCGSelective(std::move(solver), std::move(solver2), verbose, eps, maxstep,
                                 printSteps);
    if (name == "LOBPCGSelectiveFast")
      return new LOBPCGSelectiveFast(std::move(solver), std::move(solver2), verbose, eps, maxstep,
                                     printSteps);
    if (name == "LOBPCGSelectiveVeryFast")
      return new LOBPCGSelectiveVeryFast(std::move(solver), std::move(solver2), verbose, eps,
                                         maxstep, printSteps);

    THROW("EigenSolver " + name + " not implemented")
  }

  std::unique_ptr<IEigenSolver> CreateUnique() { return std::unique_ptr<IEigenSolver>(Create()); }
};

#endif // EIGENSOLVERCREATOR_HPP
