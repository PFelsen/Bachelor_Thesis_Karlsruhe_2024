#ifndef ELLIPTICPDESOLVER_HPP
#define ELLIPTICPDESOLVER_HPP

#include "EllipticProblems.hpp"
#include "IEllipticAssemble.hpp"
#include "PDESolver.hpp"
#include "Newton.hpp"

class EllipticPDESolver : public PDESolver<IEllipticProblem> {
protected:
  std::shared_ptr<Newton> newton;

  std::shared_ptr<IEllipticAssemble> assemble;

  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<IEllipticProblem> problem) override;

public:
  explicit EllipticPDESolver(const PDESolverConfig &conf);

  std::string Name() const override { return "EllipticPDESolver"; }

  void SetNormalFlux(const Vector &u, Vector &flux);

  // From PDEMain
  ValueMap ComputeValues(const Vector &u) const override;

  void PrintValues(const Solution &solution) override;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};

#endif //ELLIPTICPDESOLVER_HPP
