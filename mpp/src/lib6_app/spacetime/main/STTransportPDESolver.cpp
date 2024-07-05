#include "STTransportPDESolver.hpp"
#include "STDGDGTransportElement.hpp"
#include "SpaceTimePlotting.hpp"


void STTransportPDESolver::run(Solution &solution) const {
  solution.vector.SetAccumulateFlag(true);
  Matrix B(solution.vector);
  Vector RHS(0.0, solution.vector);
  assemble->System(B, RHS);
  (*solver)(B, true);
  RHS.SetAccumulateFlag(false);
  solution.vector = (*solver) * RHS;
  solution.converged = true;
}

void STTransportPDESolver::computeValues(Solution &solution) const {
  solution.values["L1"] = assemble->L1Norm(solution.vector);
  solution.values["L2"] = assemble->L2Norm(solution.vector);
  solution.values["Energy"] = assemble->EnergyNorm(solution.vector);
  solution.values["L2SpaceNormAtEndTime"] = assemble->L2SpaceNormAtEndTime(solution.vector);
}

void STTransportPDESolver::plotVtu(Solution &solution) const { }


void STTransportPDESolver::createAssemble(std::shared_ptr<ITransportProblem> problem) {
  problem->CreateMeshes(MeshesCreator().WithPLevel(conf.pLevel).WithLevel(conf.level));
  DegreePair degree(conf.degree, conf.degree);
  assemble = std::make_shared<STDGTransportAssemble>(degree, problem);
}
