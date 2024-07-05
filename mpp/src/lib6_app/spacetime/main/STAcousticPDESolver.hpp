#ifndef STACOUSTICPDESOLVER_HPP
#define STACOUSTICPDESOLVER_HPP

#include  "AcousticProblems.hpp"
#include "STDGViscoAcousticAssemble.hpp"
#include "SpaceTimePreconditioner.hpp"
#include "PDESolver.hpp"


class STAcousticPDESolver : public PDESolver<AcousticProblem> {
private:
  std::unique_ptr<LinearSolver> solver;

  std::shared_ptr<STDGViscoAcousticAssemble> assemble;

protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<AcousticProblem> problem) override;

public:
  explicit STAcousticPDESolver(const PDESolverConfig &conf) : GenericPDESolver(conf) {
    solver = GetLinearSolverUnique(createPreconditioner());
  }

  std::string Name() const override { return "STViscoAcousticPDESolver"; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }

private:
  void plotParams(Solution &solution) const;

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


};

#endif //STACOUSTICPDESOLVER_HPP
