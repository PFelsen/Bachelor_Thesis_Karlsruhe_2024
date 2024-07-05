#include "STAcousticPDESolver.hpp"
#include "STDGDGViscoAcousticElement.hpp"

void STAcousticPDESolver::run(Solution &solution) const {
  solution.vector.SetAccumulateFlag(true);
  Matrix B(solution.vector);
  Vector RHS(0.0, solution.vector);
  assemble->System(B, RHS);
  (*solver)(B, true);
  RHS.SetAccumulateFlag(false);
  solution.vector = (*solver) * RHS;
  solution.converged = true;
}

void STAcousticPDESolver::computeValues(Solution &solution) const {
  solution.values["L1"] = assemble->L1Norm(solution.vector);
  solution.values["L2"] = assemble->L2Norm(solution.vector);
  solution.values["Energy"] = assemble->EnergyNorm(solution.vector);
  solution.values["L2SpaceNormAtEndTime"] = assemble->L2SpaceNormAtEndTime(solution.vector);
}

void STAcousticPDESolver::plotVtu(Solution &solution) const {
//  bool parallel = true;
//  Config::Get("ParallelEstimator", parallel);
//  if (parallel) return;
//
//  if (solution.Id().coarse) return;
//  Plotting::Instance().Clear();
//  if (typeid(solution.u.GetDisc()) != typeid(STDiscretization_DGDG)) return;
//  plotParams(solution);
//
//  STDiscretizationT_DGDG vtuDisc(solution.u.GetDisc().GetMeshes(), 0, 0, 3);
//  SampleSolution solForVtu(vtuDisc, solution.Id());
//
//  for (cell c = solution.u.cells(); c != solution.u.cells_end(); ++c) {
//    SpaceTimeViscoAcousticDGTElement elem(solution.u, c, 0);
//    row r = solution.u.find_row(c());
//    VectorField velocity = elem.VelocityGlobal(c(), solution.u);
//    solForVtu.u(r, 0) = velocity[0];
//    solForVtu.u(r, 1) = velocity[1];
//    solForVtu.u(r, 2) = elem.PressureGlobal(c(), solution.u);
//  }
//
//  mpp::plot(solForVtu.IdString()).AddData("V_X", solForVtu.u, 0);
//  mpp::plot(solForVtu.IdString()).AddData("V_Y", solForVtu.u, 1);
//  mpp::plot(solForVtu.IdString()).AddData("P_0", solForVtu.u, 2);
//  mpp::plot(solForVtu.IdString()).PlotFile();
}

void STAcousticPDESolver::plotParams(Solution &solution) const {
//  STDiscretizationT_DGDG paramDisc(solution.u.GetDisc().GetMeshes(), 0, 0, 2);
//  Solution rhoAndKappa(paramDisc, solution.Id(), "Params");
//  for (row r = rhoAndKappa.u.rows(); r != rhoAndKappa.u.rows_end(); ++r) {
//    rhoAndKappa.u(r, 0) = assemble->GetProblem().Rho(r());
//    rhoAndKappa.u(r, 1) = assemble->GetProblem().Kappa(r());
//  }
//  mpp::plot(rhoAndKappa.IdString()).AddData("Rho", rhoAndKappa.u, 0);
//  mpp::plot(rhoAndKappa.IdString()).AddData("Kappa", rhoAndKappa.u, 1);
//  mpp::plot(rhoAndKappa.IdString()).PlotFile();
}

void STAcousticPDESolver::createAssemble(std::shared_ptr<AcousticProblem> problem) {
    problem->CreateMeshes(MeshesCreator().WithPLevel(conf.pLevel).WithLevel(conf.level));
    DegreePair degree((short)conf.degree, (short)conf.timeDegree);
    assemble = std::make_shared<STDGViscoAcousticAssemble>(
        problem->GetMeshes(), degree, problem
    );
}