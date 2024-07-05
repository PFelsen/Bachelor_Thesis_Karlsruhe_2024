#include "AcousticPDESolver.hpp"
#include "LinearSolver.hpp"

AcousticPDESolver::AcousticPDESolver(const PDESolverConfig &conf) : GenericPDESolver(conf) {
  timeIntegrator = TimeIntegratorCreator().
      WithLinearSolver(GetLinearSolver(GetPC())).
      WithLinearSolver(GetLinearSolver(GetPC())).
      WithLinearSolver(GetLinearSolver(GetPC())).
      WithLinearSolver(GetLinearSolver(GetPC())).
      WithRkOrder(conf.rkorder).
      CreateUnique();
}


#include "AcousticProblems.hpp"

void AcousticPDESolver::run(Solution &solution) const {
//  assemble->SetTimeSeries(solution.vector);
  solution.converged = timeIntegrator->Method(assemble.get(), solution.vector);
}

void AcousticPDESolver::computeValues(Solution &solution) const {
  solution.values["Energy"] = assemble->Energy(solution.vector);
  solution.values["L1"] = assemble->ComputeNorms(solution.vector).l1;
  solution.values["L2"] = assemble->ComputeNorms(solution.vector).l2;
  solution.values["MinP"] = assemble->MinMaxPressure(solution.vector).first;
  solution.values["MaxP"] = assemble->MinMaxPressure(solution.vector).second;
  solution.values = ComputeValues(solution.vector);
}

void AcousticPDESolver::plotVtu(Solution &solution) const {}

void AcousticPDESolver::createAssemble(std::shared_ptr<AcousticProblem> problem) {
  assemble = CreateAcousticAssembleUnique(*problem, conf);
}

std::map<std::string, double> AcousticPDESolver::ComputeValues(const Vector &u) const {
  std::map<std::string, double> values{};
  auto norms = assemble->ComputeNorms(u);

  values["L1"] = norms.l1;
  values["L2"] = norms.l2;
  values["Energy"] = norms.energy;
  Vector tmp(u);
  tmp.MakeAdditive();
  values["NormCoefficients"] = tmp.norm();

  if (!assemble->GetProblem().HasExactSolution()) return values;

  auto errors = assemble->ComputeErrors(u);

  values["L1Error"] = errors.l1;
  values["L2Error"] = errors.l2;
  values["MaxError"] = errors.inf;

  return values;
}