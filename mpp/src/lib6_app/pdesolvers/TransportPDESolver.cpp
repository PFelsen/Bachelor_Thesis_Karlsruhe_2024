#include "TransportPDESolver.hpp"
#include "GMRES.hpp"

TransportPDESolver::TransportPDESolver(const PDESolverConfig &conf) :
    timeIntegrator(TimeIntegratorCreator().
        WithLinearSolver(new GMRES(GetPC("PointBlockJacobi"))).
        WithLinearSolver(new GMRES(GetPC("PointBlockJacobi"))).
        WithLinearSolver(new GMRES(GetPC("PointBlockJacobi"))).
        WithLinearSolver(new GMRES(GetPC("PointBlockJacobi"))).
        WithConfigTypeEntry().WithRkOrder(conf.rkorder).
        CreateUnique()),
    GenericPDESolver(conf) {}


void TransportPDESolver::run(Solution &solution) const {
  TimeSeries ts;
  assemble->ResetTime(ts.FirstTStep(), ts.LastTStep(), ts.StepSize());
  solution.converged = timeIntegrator->Method(assemble.get(), solution.vector);
}

void TransportPDESolver::computeValues(Solution &solution) const {
    solution.values = ComputeValues(solution.vector);
}

void TransportPDESolver::plotVtu(Solution &solution) const {}

void TransportPDESolver::createAssemble(std::shared_ptr<ITransportProblem> problem) {
  assemble = CreateTransportAssembleUnique(*problem, conf);
}

std::map<std::string, double> TransportPDESolver::ComputeValues(const Vector &u)const  {
  std::map<std::string, double> values{};
  RatePair fluxPair = assemble->InflowOutflow();
  values["Energy"] = assemble->Energy(u);
  values["Mass"] = assemble->Mass(u);
  values["Inflow"] = fluxPair.first;
  values["Outflow"] = fluxPair.second;
  values["Flux"] = assemble->MaxFlux(u);

  if (!assemble->GetProblem().HasExactSolution()) return values;
  values["Error"] = assemble->Error(u);
  return values;
}

void TransportPDESolver::PrintValues(const Solution &solution) {
  if (verbose == 0) return;
  auto values = solution.values;
  mout << endl;
  vector<PrintInfoEntry<double>> entries{PrintInfoEntry("energy", values["Energy"]),
                                         PrintInfoEntry("total inflow", values["Inflow"]),
                                         PrintInfoEntry("total outflow", values["Outflow"]),
                                         PrintInfoEntry("mass", values["Mass"]),
                                         PrintInfoEntry("maximal flux", values["Flux"])};
  if (assemble->GetProblem().HasExactSolution()) {
    entries.emplace_back("Error", values["Error"]);
  }
  mout.PrintInfo("Transport", verbose, entries);
}

