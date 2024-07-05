#include "ReactionPDESolver.hpp"
#include "GMRES.hpp"
#include "TimeIntegrator.hpp"
#include "Newton.hpp"

void ReactionPDESolver::run(Solution &solution) const {
  //  assemble->SetTimeSeries(solution.vector);

  if (typeid(*assemble) == typeid(PGReactionAssemble)) {
    auto timeIntegrator = TimeIntegratorCreator(IMPLICIT_EULER).
        WithNonLinearSolver(new Newton(std::make_unique<GMRES>(GetPC("SuperLU")))).
        CreateNonLinearTimeIntegratorUnique();
    solution.converged = timeIntegrator->Method(*assemble, solution.vector);
  }


  if (typeid(*assemble) == typeid(DGReactionAssemble)) {
    auto timeIntegrator = TimeIntegratorCreator(IMPLICIT_EULER).
        WithNonLinearSolver(new Newton(std::make_unique<GMRES>(GetPC(
        "PointBlockJacobi")))).
        CreateNonLinearTimeIntegratorUnique();
    solution.converged = timeIntegrator->Method(*assemble, solution.vector);
  }
}

void ReactionPDESolver::computeValues(Solution &solution) const {
  solution.values["Energy"] = assemble->Energy(solution.vector);
  solution.values["Mass"] = assemble->Mass(solution.vector);
  solution.values["Inflow"] = assemble->InflowOutflow(solution.vector).first;
  solution.values["Outflow"] = assemble->InflowOutflow(solution.vector).second;
  solution.values = ComputeValues(solution.vector);
}

void ReactionPDESolver::plotVtu(Solution &solution) const {}

void ReactionPDESolver::createAssemble(std::shared_ptr<IReactionProblem> problem) {
  if (conf.modelName == "PGReaction")
    assemble = std::make_unique<PGReactionAssemble>(
        *problem, conf.degree
    );

  if (conf.modelName == "DGReaction")
    assemble = std::make_unique<DGReactionAssemble>(
        *problem, conf.degree
    );
}

std::map<std::string, double> ReactionPDESolver::ComputeValues(const Vector &u) const {
  std::map<std::string, double> values{};

  RatePair fluxPair = assemble->InflowOutflow(u);
  values["Energy"] = assemble->Energy(u);
  values["Mass"] = assemble->Mass(u);
  values["Inflow"] = fluxPair.first;
  values["Outflow"] = fluxPair.second;

  if (!assemble->GetProblem().HasExactSolution()) return values;

  return values;
}