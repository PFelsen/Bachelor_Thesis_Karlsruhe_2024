#include "ParallelSolver.hpp"
#include "Parallel.hpp"

using namespace std;

ParallelSolver::ParallelSolver() :
    min_matrix_size(20), maxP(0), PS_cd(true), steps(nullptr), PSM(0) {
  Config::Get("PSMsize", min_matrix_size);
  Config::Get("PS_maxP", maxP);
  Config::Get("PS_checkdiagonal", PS_cd);
}

void ParallelSolver::Construct(const Matrix &A) {
  Matrix AA(A);
  AA.EliminateDirichlet();
  SparseMatrix S(AA);

  if (PS_cd) S.CheckDiagonal();

  Date Start;
  steps = new ParallelSolverAllSteps(A);

  int size = steps->size();
  PSM.resize(size);
  PSM[0] = new ParallelSolverMatrix_Sparse(steps->get_step(0));
  for (int i = 1; i < size; ++i)
    PSM[i] = new ParallelSolverMatrix(steps->get_step(i), min_matrix_size, maxP);

  PSM[0]->Set(S);
  for (int i = 0; i < size - 1; ++i) {
    PSM[i]->makeLU();
    PSM[i]->SetNext_Matrix(*PSM[i + 1]);
  }
  PSM[size - 1]->makeLU();

  vout(1) << "ParallelMatrixSolver: total time (N = " << A.pSize()
          << " unknowns):" << Date() - Start << endl;
}

void ParallelSolver::Destruct() {
  for (int i = PSM.size() - 1; i >= 0; --i) {
    delete PSM[i];
  }
  PSM.clear();
  if (steps) delete steps;
  steps = nullptr;
}

void ParallelSolver::multiply(Vector &u, const Vector &b) const {
  int size = steps->size();

  for (int i = 0; i < size; ++i)
    PSM[i]->Create_rhs();
  PSM[0]->Set_rhs(b());

  for (int i = 0; i < size - 1; ++i) {
    PSM[i]->SolveL();
    PSM[i]->SetNext_rhs_LEFT(*PSM[i + 1]);
  }
  PSM[size - 1]->SolveL();

  for (int i = size - 1; i >= 1; --i) {
    PSM[i]->SolveU();
    PSM[i]->SetNext_rhs_RIGHT(*PSM[i - 1]);
  }
  PSM[0]->SolveU();
  PSM[0]->Write_rhs(u());
  u.SetAccumulateFlag(true); // necessary since u is set component wise
}

void ParallelSolver::multiply(Vectors &U, const Vectors &B) const {
  int size = steps->size();

  int nrhs = U.size();
  for (int i = 0; i < size; ++i)
    PSM[i]->Create_rhs(nrhs);

  for (int i = 0; i < nrhs; ++i)
    PSM[0]->Set_rhs(B[i](), i);

  for (int i = 0; i < size - 1; ++i) {
    PSM[i]->SolveL();
    PSM[i]->SetNext_rhs_LEFT(*PSM[i + 1]);
  }
  PSM[size - 1]->SolveL();

  for (int i = size - 1; i >= 1; --i) {
    PSM[i]->SolveU();
    PSM[i]->SetNext_rhs_RIGHT(*PSM[i - 1]);
  }
  PSM[0]->SolveU();

  for (int i = 0; i < nrhs; ++i)
    PSM[0]->Write_rhs(U[i](), i);
}
