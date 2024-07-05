#include "PDEModel.hpp"
#include "LagrangeEllipticAssemble.hpp"

EllipticInverseProblem &EllipticPDEModel::GetProblem() { return *problem; }

EllipticPDEModel::EllipticPDEModel(PDESolverConfig conf) :
    EllipticPDESolver(conf), observation_operator(conf.observations) {
  problem = CreateEllipticInverseShared(conf.problemName, conf.proposal_config);
  problem->CreateMeshes(MeshesCreator().WithPLevel(1).WithLevel(2));
//  problem->GetInput().proposal0D = 1.0;
}

void EllipticPDEModel::Run(RVector &solution) { run(solution); }

bool EllipticPDEModel::run(RVector &observation) {
  mout.StartBlock("PDEModel");
  assemble = CreateEllipticAssembleShared(*problem, conf);
  Vector solution(0.0, assemble->GetSharedDisc());
  mout << "Start solving on l=" << solution.Level().space << endl;
  solution.SetAccumulateFlag(true);
  newton->operator()(*assemble, solution);
  PPM->Barrier(0);
  // Is this really the observation or rather the current guess evaluated at the observation points
  observation = observation_operator.GetObservation(solution);
  // pout << "End evaluate observation" << endl;
  mout.EndBlock();
  return true;
}

TestModel::TestModel(const PDESolverConfig &conf) : EllipticPDEModel(conf) {}

void TestModel::Run(RVector &solution) {
  solution = problem->GetInput().GetRandomField().GetParameterVector();
}

EllipticPDEModel *CreatePDEModel(PDESolverConfig conf) {
  if (conf.modelName.find("Elliptic")) { return new EllipticPDEModel(conf); }
  if (conf.modelName.find("Test")) { return new TestModel(conf); }
  throw std::invalid_argument("PDE-model type " + conf.modelName + " cannot be found.");
}

std::shared_ptr<EllipticPDEModel> CreatePDEModelShared(PDESolverConfig conf) {
  return std::shared_ptr<EllipticPDEModel>(CreatePDEModel(std::move(conf)));
}
