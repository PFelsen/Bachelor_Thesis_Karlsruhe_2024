#ifndef ACOUSTICPDESOLVER_HPP
#define ACOUSTICPDESOLVER_HPP

#include "AcousticProblems.hpp"
#include "IAcousticAssemble.hpp"
#include "PDESolver.hpp"
#include "TimeIntegrator.hpp"

class AcousticPDESolver : public PDESolver<AcousticProblem> {
private:
  std::shared_ptr<IAcousticAssemble> assemble;

  std::unique_ptr<TimeIntegrator> timeIntegrator;

protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<AcousticProblem> problem) override;

public:
  explicit AcousticPDESolver(const PDESolverConfig &conf);

  std::string Name() const override { return "AcousticPDESolver"; }

  ValueMap ComputeValues(const Vector &u) const override;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};

#endif //ACOUSTICPDESOLVER_HPP
