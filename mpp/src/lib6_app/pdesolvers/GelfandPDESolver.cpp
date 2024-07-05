#include "GelfandPDESolver.hpp"

GelfandPDESolver::GelfandPDESolver(const PDESolverConfig &conf) :
    newton(CreateNewton(conf.linearSolver, conf.preconditioner)), GenericPDESolver(conf) {}

void GelfandPDESolver::run(Solution &solution) const {
  solution.vector.SetAccumulateFlag(true);
  assemble->Initialize(solution.vector);
  newton->operator()(*assemble, solution.vector);
  solution.converged = newton->converged();
}

void GelfandPDESolver::computeValues(Solution &solution) const {
  solution.values = ComputeValues(solution.vector);
}

void GelfandPDESolver::plotVtu(Solution &solution) const {
  if (!plotting) return;
  mpp::plot("u") << solution.vector << mpp::endp;
}

void GelfandPDESolver::createAssemble(std::shared_ptr<GelfandProblem> problem) {
  try {
    assemble =
        CreateGelfandAssembleUnique(*std::dynamic_pointer_cast<GelfandProblem>(problem), conf);
  } catch (const std::bad_cast &ex) { THROW("PDESolver can't solve problems of given type") }
}

std::map<std::string, double> GelfandPDESolver::ComputeValues(const Vector &u) const {
  std::map<std::string, double> values{};

  values["DoFCount"] = u.size();
  values["LinearSteps"] = newton->GetLinearSolver().GetIteration().Steps();

  values["Energy"] = assemble->Energy(u);

  return values;
}

void GelfandPDESolver::PrintValues(const Solution &solution) {
  if (verbose == 0) return;
  else mout << "max(u) = " << solution.vector.Max() << endl;
}