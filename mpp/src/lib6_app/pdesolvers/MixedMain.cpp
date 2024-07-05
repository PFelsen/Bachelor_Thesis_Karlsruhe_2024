#include "HybridEllipticAssemble.hpp"
#include "MixedEllipticAssemble.hpp"
#include "MeshesCreator.hpp"
#include "HybridMain.hpp"
#include "MixedMain.hpp"
#include "Newton.hpp"
#include "LagrangeDiscretization.hpp"

template<typename Assemble>
void plotVtu(std::shared_ptr<Assemble> assemble, std::shared_ptr<IProblem> problem, const Vector &u) {
  auto  fluxDisc = std::make_shared<const LagrangeDiscretization>(problem->GetMeshes(), 0, 3);
  auto  uDisc = std::make_shared<const LagrangeDiscretization>(problem->GetMeshes(), 0, 1);
  Vector flux(0.0, fluxDisc);
  Vector p(0.0, uDisc);
  assemble->SetPressureFlux(u, p, flux);

  mpp::plot("Flux") << flux << mpp::endp;
  mpp::plot("P") << p << mpp::endp;
}

std::shared_ptr<IEllipticProblem> setupProblem(const std::string &problemName) {
  auto problem = CreateEllipticProblemShared(problemName);
  return problem;
}

MixedMain::MixedMain(const std::string &problemName) {
  problem = setupProblem(problemName);
  assemble = std::make_shared<MixedEllipticAssemble>(*problem);
  const Meshes &M = assemble->GetSharedDisc()->GetMeshes();

  for (LevelPair level = M.PLevel();
        level != M.FineLevel().NextInSpace();
        level = level.NextInSpace()) {
    solutionOnLevel[level.space] = std::make_shared<Vector>(0.0, assemble->GetSharedDisc(), level);
  }
  
  Config::Get("MixedVerbose", verbose);
  Solve();
}

void MixedMain::PlotVtu(const Vector &u) const {
  plotVtu(assemble, problem, u);
}

IProblem& MixedMain::GetFluxProblem() {
  return *problem;
}

VectorField MixedMain::EvaluateCellFlux(const LevelPair& level, const Cell &c) const {
  return assemble->EvaluateCellFlux(*solutionOnLevel.at(level.space), c.SpaceCell());
}

double MixedMain::EvaluateNormalFlux(const LevelPair& level, const Cell &c, int face) const {
  return assemble->EvaluateNormalFlux(*solutionOnLevel.at(level.space), c.SpaceCell(), face);
}

void MixedMain::Solve() {
  assemble->GetSharedDisc()->GetMeshes().PrintInfo();
  mout.StartBlock("Mixed Flux Problem");

  for (auto &[spacelevel, solution] : solutionOnLevel) {
    vout(1) << "Solve on level: " << spacelevel << endl;
    NewtonMethod(*assemble, *solution, true);
  }
  PlotVtu(*solutionOnLevel.at(assemble->GetSharedDisc()->GetMeshes().FineLevel().space));

  
  mout.EndBlock(verbose < 1);
  vout(1) << endl;
}

HybridMain::HybridMain(const std::string &problemName) {
  problem = setupProblem(problemName);
  assemble = std::make_shared<HybridEllipticAssemble>(*problem);
  solution = std::make_shared<Vector>(assemble->GetSharedDisc());
  normalFlux = std::make_shared<Vector>(assemble->GetSharedDisc());

  Config::Get("HybridVerbose", verbose);
  Solve();
}

void HybridMain::PlotVtu(const Vector &u) const {
  plotVtu(assemble, problem, u);
}

const IProblem& HybridMain::GetFluxProblem() const {
    return *problem;
}

VectorField HybridMain::EvaluateCellFlux(const Cell &c) const {
  return assemble->EvaluateCellFlux(*normalFlux, c.SpaceCell());
}

double HybridMain::EvaluateNormalFlux(const Cell &c, int face) const {
  return assemble->EvaluateNormalFlux(*normalFlux, c.SpaceCell(), face);
}

void HybridMain::Solve() {
  assemble->GetSharedDisc()->GetMeshes().PrintInfo();
  mout.StartBlock("Hybrid Flux Problem");
  vout(1) << "Solve" << endl;
  *solution = 0.0;
  NewtonMethod(*assemble, *solution, true);
  assemble->SetNormalFlux(*solution, *normalFlux);

  PlotVtu(*solution);

  mout.EndBlock(verbose < 1);
  vout(1) << endl;
}
