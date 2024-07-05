#include "EllipticPDESolver.hpp"
#include "HybridEllipticAssemble.hpp"
#include "LagrangeEllipticAssemble.hpp"
#include "MixedEllipticAssemble.hpp"

EllipticPDESolver::EllipticPDESolver(const PDESolverConfig &conf) :
    newton(CreateNewton(conf.linearSolver, conf.preconditioner)), GenericPDESolver(conf) {}

void EllipticPDESolver::run(Solution &solution) const {
  solution.vector.SetAccumulateFlag(true);
  newton->operator()(*assemble, solution.vector);
  solution.converged = newton->converged();
}

void EllipticPDESolver::computeValues(Solution &solution) const {
  solution.values = ComputeValues(solution.vector);
}

void EllipticPDESolver::plotVtu(Solution &solution) const {
  if (!plotting) return;

  auto &meshes = assemble->GetProblem().GetMeshes();

  std::shared_ptr<IDiscretization> pressureDisc;
  if (typeid(*assemble) == typeid(HybridEllipticAssemble)
      || typeid(*assemble) == typeid(MixedEllipticAssemble)) {
    pressureDisc = std::make_shared<LagrangeDiscretization>(meshes, 0, 1);
  } else {
    pressureDisc = std::make_shared<LagrangeDiscretization>(meshes, 1, 1);
  }

  auto fluxDisc = std::make_shared<const LagrangeDiscretization>(meshes, 0, SpaceDimension);

  Vector pressure(0.0, pressureDisc);
  Vector flux(0.0, fluxDisc, solution.vector.Level());
  if (typeid(*assemble) == typeid(HybridEllipticAssemble)) {
    auto hybridAssemble = std::dynamic_pointer_cast<HybridEllipticAssemble>(assemble);

    hybridAssemble->SetPressureFlux(solution.vector, pressure, flux);
  } else {
    assemble->SetPressure(solution.vector, pressure);
    assemble->SetFlux(solution.vector, flux);
  }

  mpp::plot("u") << pressure << mpp::endp;
  mpp::plot("flux") << flux << mpp::endp;

  auto kappaDisc = std::make_shared<const LagrangeDiscretization>(meshes, 0, 1);
  Vector kappa(0.0, kappaDisc, solution.vector.Level());
  assemble->GetProblem().Permeability(kappa);
  mpp::plot("kappa") << kappa << mpp::endp;

  if (assemble->GetProblem().HasExactSolution()) {
    Vector u_ex(0.0, solution.vector);
    assemble->SetExactSolution(u_ex);
    Vector error = solution.vector - u_ex;

    if (typeid(*assemble) == typeid(HybridEllipticAssemble)) {
      Vector u_ex_plottable(0.0, pressureDisc);
      Vector u_error_plottable(0.0, pressureDisc);
      auto hybridAssemble = std::dynamic_pointer_cast<HybridEllipticAssemble>(assemble);
      hybridAssemble->SetPressureFlux(u_ex, u_ex_plottable, flux);
      mpp::plot("u_ex") << u_ex_plottable << mpp::endp;
      hybridAssemble->SetPressureFlux(error, u_error_plottable, flux);
      mpp::plot("u_error") << u_error_plottable << mpp::endp;
    } else {
      mpp::plot("u_ex") << u_ex << mpp::endp;
      mpp::plot("u_error") << error << mpp::endp;
    }
  }

  //  if (!plotting) return;
  //  if (typeid(*assemble).name() == typeid(HybridEllipticAssemble).name()) return;
  //  if (typeid(*assemble).name() == typeid(MixedEllipticAssemble).name()) return;
  //
  //  LagrangeDiscretization pressureDisc(assemble->GetProblem().GetMeshes(), 0, 1);
  //  Solution pressure(pressureDisc);
  //
  //  assemble->SetPressure(solution.vector, pressure.vector);
  //
  //  mpp::PlotData("Pressure", "Pressure" + assemble->GetProblem().Name(), solution.vector);
  //
  //  LagrangeDiscretization fluxDisc(GetDisc().GetMeshes(), 0, SpaceDimension);
  //  Solution flux(fluxDisc, solution.vector.Level());
  //  assemble->SetFlux(solution.vector, flux.vector);
  //  mpp::PlotData("Flux", "Flux" + assemble->GetProblem().Name(), flux.vector);
  //
  //  LagrangeDiscretization kappaDisc(GetDisc().GetMeshes(), 0, 1);
  //  Solution kappa(kappaDisc, solution.vector.Level());
  //  assemble->IEllipticAssemble::GetProblem().Permeability(kappa.vector);
  //  mpp::PlotData("Kappa", "Kappa" + assemble->GetProblem().Name(), kappa.vector);
  //
  //  if (assemble->GetProblem().HasExactSolution()) {
  //    Solution exactSolution(solution.vector);
  //    assemble->SetExactSolution(exactSolution.vector);
  //    mpp::plot("u_ex") << exactSolution.vector << mpp::endp;
  //
  //    Vector error = solution.vector - exactSolution.vector;
  //    mpp::plot("u_error") << error << mpp::endp;
  //  }
}

void EllipticPDESolver::SetNormalFlux(const Vector &u, Vector &flux) {
  if (abs(assemble->FluxError(u)) < 1e-10) assemble->SetNormalFlux(u, flux);
  else Warning("Flux Error is not Zero but used as Normal Flux.")
}

void EllipticPDESolver::createAssemble(std::shared_ptr<IEllipticProblem> problem) {
  try {
    assemble =
        CreateEllipticAssembleUnique(*std::dynamic_pointer_cast<IEllipticProblem>(problem), conf);
  } catch (const std::bad_cast &ex) { THROW("PDESolver can't solve problems of given type") }
}

std::map<std::string, double> EllipticPDESolver::ComputeValues(const Vector &u) const {
  std::map<std::string, double> values{};

  values["DoFCount"] = u.size();
  values["LinearSteps"] = newton->GetLinearSolver().GetIteration().Steps();
  FluxPair fluxPair = assemble->InflowOutflow(u);
  FluxPair preFluxPair = assemble->PrescribedInflowOutflow(u);

  values["L2"] = assemble->L2(u);
  values["H1"] = assemble->H1(u);
  values["Energy"] = assemble->Energy(u);
  values["Inflow"] = fluxPair.first;
  values["Outflow"] = fluxPair.second;
  values["PreInflow"] = preFluxPair.first;
  values["PreOutflow"] = preFluxPair.second;
  values["Goal"] = assemble->GoalFunctional(u);
  values["FluxLoss"] = fluxPair.first + fluxPair.second;
  values["FluxError"] =
      assemble->FluxError(u); // Checks neumann bc, not dependent of HasExactSolution

  if (conf.dualPrimal) values["DualPrimal"] = assemble->DualPrimalError(u);

  if (assemble->GetProblem().Name() == "Rock") {
    FluxPair outflowLeftRight = assemble->OutflowLeftRight(u);
    values["OutflowLeft"] = outflowLeftRight.first;
    values["OutflowRight"] = outflowLeftRight.second;
  }

  if (!assemble->GetProblem().HasExactSolution()) return values;

  values["L2Error"] = assemble->L2Error(u);
  values["MaxError"] = assemble->MaxError(u);
  values["FaceError"] = assemble->FaceError(u);
  values["EnergyError"] = assemble->EnergyError(u);
  values["L2CellAvgError"] = assemble->L2CellAvgError(u);

  return values;
}

void EllipticPDESolver::PrintValues(const Solution &solution) {
  if (verbose == 0) return;

  auto values = solution.values.empty() ? ComputeValues(solution.vector) : solution.values;

  mout << endl;


  mout.PrintInfo("Mixed", verbose, PrintInfoEntry("DoFCount", values["DoFCount"]));


  mout.PrintInfo("Mixed", verbose, PrintInfoEntry("DoFCount", values["DoFCount"]),
                 PrintInfoEntry("LinearSteps", values["LinearSteps"]),
                 PrintInfoEntry("Prescribed Inflow", values["PreInflow"]),
                 PrintInfoEntry("Prescribed Outflow", values["PreOutflow"]),
                 PrintInfoEntry("Calculated Inflow", values["Inflow"]),
                 PrintInfoEntry("Calculated Outflow", values["Outflow"]),
                 PrintInfoEntry("Flux Loss", values["FluxLoss"]),
                 PrintInfoEntry("H1 Norm", values["H1"]), PrintInfoEntry("L2 Norm", values["L2"]),
                 PrintInfoEntry("Goal Functional", values["Goal"]));

  if (assemble->GetProblem().Name() == "Rock") {
    mout.PrintInfo("Rock Problem", verbose, PrintInfoEntry("Left Outflow", values["OutflowLeft"]),
                   PrintInfoEntry("Right Outflow", values["OutflowRight"]));
  }

  if (conf.dualPrimal)
    mout.PrintInfo("Dual Primal Error", verbose,
                   PrintInfoEntry("Dual Primal Error", values["DualPrimal"]));

  if (assemble->GetProblem().HasExactSolution()) {
    mout.PrintInfo("Exact Solution", verbose, PrintInfoEntry("Supremum Error", values["MaxError"]),
                   PrintInfoEntry("Energy Error", values["EnergyError"]),
                   PrintInfoEntry("L2 Error", values["L2Error"]),
                   PrintInfoEntry("L2 Error Average", values["L2CellAvgError"]),
                   PrintInfoEntry("Face Error", values["FaceError"]),
                   PrintInfoEntry("Flux Error", values["FluxError"]));
  }
}
