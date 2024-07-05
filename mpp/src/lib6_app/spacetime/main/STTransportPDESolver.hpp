#ifndef MPP_STTRANSPORTPDESOLVER_H
#define MPP_STTRANSPORTPDESOLVER_H

#include "STTransportProblems.hpp"
#include "STDGTransportAssemble.hpp"
#include "SpaceTimePreconditioner.hpp"
#include "PDESolver.hpp"


class STTransportPDESolver : public PDESolver<ITransportProblem> {
private:
  std::unique_ptr<LinearSolver> solver;

  std::shared_ptr<STDGTransportAssemble> assemble;

  std::unique_ptr<Preconditioner> createPreconditioner() {
    std::string pathChoice;
    Config::Get("PathChoice", pathChoice);
    if (pathChoice == "space_time") {
      return std::make_unique<STMultiGridPC>(
          std::make_unique<SpaceThenTimePathStrategy>(), *assemble
      );
    }
    else if (pathChoice == "time_space") {
      return std::make_unique<STMultiGridPC>(
          std::make_unique<TimeThenSpacePathStrategy>(), *assemble
      );
    }
    return std::unique_ptr<Preconditioner>(GetPC("PointBlockGaussSeidel"));
  }

protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<ITransportProblem> problem) override;

public:
  explicit STTransportPDESolver(const PDESolverConfig &conf) : GenericPDESolver(conf) {
    solver = GetLinearSolverUnique(createPreconditioner());
  }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }

  std::string Name() const override { return "STTransportPDESolver"; }

};

#endif  // MPP_STTRANSPORTPDESOLVER_H
