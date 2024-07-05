#ifndef TRANSPORTPDESOLVER_HPP
#define TRANSPORTPDESOLVER_HPP

#include "TransportProblems.hpp"
#include "ITransportAssemble.hpp"
#include "PDESolver.hpp"
#include "TimeIntegrator.hpp"

class TransportPDESolver : public PDESolver<ITransportProblem> {
private:
  std::shared_ptr<ITransportAssemble> assemble;

  std::unique_ptr<TimeIntegrator> timeIntegrator;

protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<ITransportProblem> problem) override;

public:
  explicit TransportPDESolver(const PDESolverConfig &conf);

  std::string Name() const override { return "TransportPDESolver"; }

  ValueMap ComputeValues(const Vector &u) const override;

  void PrintValues(const Solution &solution) override;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};

#endif //TRANSPORTPDESOLVER_HPP
