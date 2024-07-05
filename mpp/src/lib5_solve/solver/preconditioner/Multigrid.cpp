#include "Multigrid.hpp"

#include "Matrix.hpp"
#include "Vector.hpp"

const Matrix *Multigrid::StepData::createMatrix(const Vector &v, const IAssemble &assemble) {
  Matrix *M = new Matrix(v);
  assemble.Jacobi(v, *M);
  return M;
}

Multigrid::StepData::StepData(const Vector &v, const IAssemble &assemble) :
    StepData(createMatrix(v, assemble), assemble) {}

Multigrid::StepData::StepData(const Matrix *a, const IAssemble &assemble) : A(a) {
  Config::Get("presmoothing", presmooth);
  Config::Get("postsmoothing", postsmooth);
  string smootherName = "Jacobi";
  Config::Get("Smoother", smootherName);
  smoother = std::unique_ptr<Preconditioner>(GetPC(smootherName));

  string cname = "V";
  Config::Get("cycle", cname);
  if (cname == "W") {
    cycle = 2;
  } else if (cname == "0") {
    cycle = 0;
  } else if (cname == "VW") {
    cycle = 3;
  }

  auto disc = A->GetSharedDisc();
  LevelPair level = A->Level();
  LevelPair level_coarse = level.CoarserInSpace();
  LevelPair fine_level = A->GetDisc().GetMeshes().FineLevel();
  int dynamicSmoothing = 0;
  Config::Get("DynamicSmoothing", dynamicSmoothing);

  smoother->Construct(*A);
  v_coarse = std::make_unique<Vector>(disc, level_coarse);
  assemble.Initialize(*v_coarse);

  transfer = GetTransfer(*v_coarse, A->GetVector());
  // transfer->Project(*v_coarse, A->GetVector());
  *v_coarse = A->GetVector() * (*transfer);
}

void Multigrid::Construct(const Matrix &) {
  THROW("Construct without Assemble not possible for MultigridPreconditioner")
}

void Multigrid::Construct(const Matrix &A, const IAssemble &assemble) {
  fineLevel = A.Level();
  baseLevel = A.GetDisc().GetMeshes().CoarseLevel();

  if (fineLevel != baseLevel) {
    if (stepData.find(fineLevel) != stepData.end() && stepData[fineLevel].A == &A) return;
    if (!stepData.empty()) { Destruct(); }
    stepData[fineLevel] = StepData(&A, assemble);

    for (LevelPair l = fineLevel; l != baseLevel.NextInSpace(); l = l.CoarserInSpace()) {
      stepData[l.CoarserInSpace()] = StepData(*stepData[l].v_coarse, assemble);
    }
    Vector &baseV = *stepData[baseLevel.NextInSpace()].v_coarse;
    baseA = std::make_unique<Matrix>(baseV);
    assemble.Jacobi(baseV, *baseA);
  } else {
    baseRHS = std::make_unique<Vector>(A.GetSharedDisc(), fineLevel);
    assemble.Initialize(*baseRHS);
    baseA = std::make_unique<Matrix>(*baseRHS);
    assemble.Jacobi(*baseRHS, *baseA);
  }
  baseSolver = std::unique_ptr<LinearSolver>(GetLinearSolverByPrefix("Base"));
  (*baseSolver)(*baseA);
  if (!Config::Exists("BaseVerbose")) { baseSolver->SetVerbose(-1); }
}

void Multigrid::Destruct() {
  for (auto &[key, data] : stepData) {
    data.smoother->Destruct();
    if (key != fineLevel) { delete data.A; }
  }
  stepData.clear();
}

void Multigrid::Cycle(Vector &u, Vector &r) const {
  if (u.Level() == baseLevel) {
    vout(2).StartBlock("BaseSolver");
    vout(2) << "Solving on level " << u.Level() << endl;
    u = (*baseSolver) * r;
    if (u.identify()) {
      u.MakeAdditive();
      u.Accumulate();
    }
    vout(2).EndBlock(true);
    return;
  }
  const StepData &data = stepData.at(u.Level());

  vout(2).StartBlock("Multigrid");
  vout(2) << data.presmooth << " presmoothing-steps on level " << u.Level() << " with ["
          << data.smoother->Name() << ", Damp: " << theta << "]" << endl;

  Vector w(u);
  for (int i = 0; i < data.presmooth; ++i) {
    w = (*data.smoother) * r;
    if (theta != 1) w *= theta;
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
          << data.smoother->Name() << ", Damp: " << theta << "]" << endl;
  for (int i = 0; i < data.postsmooth; ++i) {
    w = r * (*data.smoother);
    if (theta != 1) w *= theta;
    r -= (*data.A) * w;
    u += w;
  }
  vout(2).EndBlock(true);
}

void Multigrid::multiply(Vector &u, const Vector &b) const {
  Vector r = b;
  u = 0;
  Cycle(u, r);
}

void Multigrid::multiply_transpose(Vector &u, Vector &v) const { multiply(u, v); }

void Multigrid::Cycle(Vectors &u, Vectors &r) const {
  if (u.Level() == baseLevel) {
    vout(2).StartBlock("BaseSolver");
    vout(2) << "Solving on level " << u.Level() << endl;
    u = (*baseSolver) * r;
    if (u.identify()) {
      u.MakeAdditive();
      u.Accumulate();
    }
    vout(2).EndBlock(true);
    return;
  }
  const StepData &data = stepData.at(u.Level());

  vout(2).StartBlock("Multigrid");
  vout(2) << data.presmooth << " presmoothing-steps on level " << u.Level() << " with ["
          << data.smoother->Name() << ", Damp: " << theta << "]" << endl;

  Vectors w(u);
  Vectors d(int(r.size()), *data.v_coarse);
  Vectors c(d);
  for (int i = 0; i < data.presmooth; ++i) {
    w = (*data.smoother) * r;
    if (theta != 1) w *= Scalar(theta);
    r -= (*data.A) * w;
    u += w;
  }


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
          << data.smoother->Name() << ", Damp: " << theta << "]" << endl;
  for (int i = 0; i < data.postsmooth; ++i) {
    w = r * (*data.smoother);
    if (theta != 1) w *= Scalar(theta);
    r -= (*data.A) * w;
    u += w;
  }
  vout(2).EndBlock(true);
}

void Multigrid::multiply(Vectors &u, const Vectors &b) const {
  mout.StartBlock(Name());
  Vectors r = b;
  u = 0;
  Cycle(u, r);
  mout.EndBlock(verbose <= 0);
}

void Multigrid::multiply_transpose(Vectors &u, Vectors &v) const { multiply(u, v); }