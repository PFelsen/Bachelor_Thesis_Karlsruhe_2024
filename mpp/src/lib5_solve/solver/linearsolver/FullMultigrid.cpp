#include "FullMultigrid.hpp"

#include "Matrix.hpp"
#include "Vector.hpp"

FullMultigrid::StepData::StepData(FullMultigrid *fm, const IAssemble &assemble,
                                  std::shared_ptr<const IDiscretization> disc, LevelPair level) :
    assemble(assemble) {
  LevelPair baseLevel = disc->GetMeshes().CoarseLevel();

  string cycleName = "V";
  Config::Get("cycle", cycleName);
  if (cycleName == "W") {
    cycle = 2;
  } else if (cycleName == "0") {
    cycle = 0;
  } else if (cycleName == "VW") {
    cycle = 3;
  }

  LevelPair level_coarse = level.CoarserInSpace();
  solution = std::make_unique<Vector>(0.0, disc, level);
  assemble.Initialize(*solution);
  rhs = std::make_unique<Vector>(0.0, *solution);
  assemble.Residual(*solution, *rhs);
  A = std::make_unique<Matrix>(*solution);
  assemble.Jacobi(*solution, *A);

  Config::Get("presmoothing", presmooth);
  Config::Get("postsmoothing", postsmooth);
  string smootherName = "Jacobi";
  Config::Get("Smoother", smootherName);
  smoother = std::unique_ptr<Preconditioner>(GetPC(smootherName));
  Config::Get("SmootherDamp", theta);
  smoother->Construct(*A);

  string solverName = "gmres";
  if (level == baseLevel) {
    Config::Get("BaseSolver", solverName);
    solver = std::unique_ptr<LinearSolver>(GetLinearSolverByPrefix("Base"));
    if (!Config::Exists("BaseVerbose")) { solver->SetVerbose(-1); }
  } else {
    Config::Get("FMLinearSolver", solverName);
    solver = GetLinearSolverUnique(solverName, std::make_unique<MultigridPreconditioner>(fm),
                                   "FMLinear");
    v_coarse = std::make_unique<Vector>(disc, level_coarse);
    assemble.Initialize(*v_coarse);
    transfer = GetTransfer(*v_coarse, *rhs);
    transfer->Project(*v_coarse, *rhs);
  }
  (*solver)(*A, assemble);
}

void FullMultigrid::MultigridPreconditioner::Cycle(Vector &u, Vector &r) const {
  vout(3) << "Calling Cycle at level " << u.Level().space << endl;
  const StepData &data = fm->stepData.at(u.Level());
  if (u.Level() == fm->baseLevel) {
    vout(2).StartBlock("BaseSolver");
    vout(2) << "Solving on level " << u.Level() << endl;
    u = (*data.solver) * r;
    if (u.identify()) {
      u.MakeAdditive();
      u.Accumulate();
    }
    vout(2).EndBlock(true);
    return;
  }


  vout(2).StartBlock("Multigrid");
  vout(2) << data.presmooth << " presmoothing-steps on level " << u.Level() << " with ["
          << data.smoother->Name() << ", Damp: " << data.theta << "]" << endl;

  Vector w(u);
  for (int i = 0; i < data.presmooth; ++i) {
    w = (*data.smoother) * r;
    if (data.theta != 1) w *= data.theta;
    r -= (*data.A) * w;
    u += w;
  }

  Vector d(*data.v_coarse);
  Vector c(d);
  for (int i = 0; i < data.cycle; ++i) {
    // ProlongateTransposed
    d = r * (*data.transfer);
    c = 0;
    Cycle(c, d);
    // Prolongate
    w = (*data.transfer) * c;
    u += w;
    r -= (*data.A) * w;
  }

  vout(2) << data.postsmooth << " postsmoothing-steps on level " << u.Level() << " with ["
          << data.smoother->Name() << ", Damp: " << data.theta << "]" << endl;
  for (int i = 0; i < data.postsmooth; ++i) {
    w = r * (*data.smoother);
    if (data.theta != 1) w *= data.theta;
    r -= (*data.A) * w;
    u += w;
  }
  vout(2).EndBlock(true);
}

void FullMultigrid::multiply(Vector &u, const Vector &b) const {
  bool reusePreviousSolution = true;
  Config::Get("ReusePreviousSolution", reusePreviousSolution);
  std::unique_ptr<Vector> previousSolution = nullptr;
  for (LevelPair level = baseLevel; level != fineLevel.NextInSpace(); level = level.NextInSpace()) {
    mout.StartBlock("FullMultigridOnLevel" + std::to_string(level.space));
    vout(3) << "Solving on level " << level.space << endl;
    const StepData &data = stepData.at(level);
    Vector &rhs = *data.rhs;
    Vector &solution = *data.solution;
    Matrix &A = *data.A;
    ITransfer &transfer = *data.transfer;
    if (reusePreviousSolution && previousSolution) {
      Vector &prevSolution = *previousSolution;
      solution = transfer * prevSolution;
      rhs -= A * solution;
    }
    solution += (*data.solver) * rhs;
    previousSolution = std::make_unique<Vector>(solution);
    mout.EndBlock();
  }
  u = *previousSolution;
}
