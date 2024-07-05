//
// Created by lstengel on 21.07.22.
//

#include "VectorValuedMain.hpp"
#include "IEllipticAssemble.hpp"


void VectorValuedMain::createAssemble(std::shared_ptr<IVectorValuedProblem> problem) {
  assemble = CreateVectorValuedAssembleShared(*problem, conf);
}

void VectorValuedMain::run(Solution &solution) const {
  solution.vector.SetAccumulateFlag(true);
  newton->operator()(*assemble, solution.vector);
  solution.converged = newton->converged();
}


void VectorValuedMain::plotVtu(Solution &solution) const { }

void VectorValuedMain::computeValues(Solution &solution) const  {
  const Vector &u = solution.vector;
  std::map<std::string, double> &values = solution.values;

  values["L2"] = assemble->L2(u);
  values["H1"] = assemble->H1(u);
  values["Energy"] = assemble->Energy(u);

  if (!assemble->GetProblem().HasExactSolution()) return;

  values["L2Error"] = assemble->L2Error(u);
  values["MaxError"] = assemble->MaxError(u);
  values["FaceError"] = assemble->FaceError(u);
  values["FluxError"] = assemble->FluxError(u);
  values["EnergyError"] = assemble->EnergyError(u);
  values["L2CellAvgError"] = assemble->L2CellAvgError(u);
}
double VectorValuedMain::EvaluateQuantity(const Vector &u, const std::string &quantity) const {
  if (quantity == "L2") { return assemble->L2Error(u); }
  else if (quantity == "Energy") { return assemble->EnergyError(u); }
  else { return infty; }
}

void VectorValuedMain::PrintValues(const Solution &solution) {
  if (verbose == 0) return;

  auto values = solution.values.empty() ? ComputeValues(solution.vector) : solution.values;

  mout << endl;


  mout.PrintInfo("Solution Measures", verbose,
                 PrintInfoEntry("H1 Norm", values["H1"]),
                 PrintInfoEntry("L2 Norm", values["L2"]));

  mout.PrintInfo("Flux", verbose,
                 PrintInfoEntry("Flux Loss", values["FluxLoss"]),
                 PrintInfoEntry("Flux Error", values["FluxError"]));

  if (assemble->GetProblem().HasExactSolution()) {
    mout.PrintInfo("Exact Solution", verbose,
                   PrintInfoEntry("Supremum Error", values["MaxError"]),
                   PrintInfoEntry("Energy Error", values["EnergyError"]),
                   PrintInfoEntry("L2 Error", values["L2Error"]),
                   PrintInfoEntry("L2 Error Average", values["L2CellAvgError"]),
                   PrintInfoEntry("Face Error", values["FaceError"]));
  }
}
