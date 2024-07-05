#include "Assemble.hpp"

double IAssemble::Energy(const Vector &u) const {
  double energy = 0;
  TRY {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      Energy(c, u, energy);
    }
  }
  CATCH("Error in Energy")
  return sqrt(PPM->SumOnCommSplit(energy, u.CommSplit()));
}

double IAssemble::Residual(const Vector &u, Vector &defect) const {
  defect = 0;
  TRY {
    for (cell ce = u.cells(); ce != u.cells_end(); ++ce)
      Residual(ce, u, defect);
  }
  CATCH("Error in Residual")
  defect.ClearDirichletValues();
  defect.Collect();
  return defect.norm();
}

void IAssemble::Jacobi(const Vector &u, Matrix &jacobi) const {
  jacobi = 0;
  TRY {
    for (cell c = jacobi.cells(); c != jacobi.cells_end(); ++c)
      Jacobi(c, u, jacobi);
  }
  CATCH("Error in Jacobi")
  jacobi.ClearDirichletValues();
}

void ILinearAssemble::AssembleSystem(Matrix &systemMatrix, Vector &rhs) const {
  rhs = 0;
  systemMatrix = 0;
  TRY {
    for (cell ce = rhs.cells(); ce != rhs.cells_end(); ++ce)
      AssembleSystem(ce, systemMatrix, rhs);
  }
  CATCH("Error in Residual")
  rhs.ClearDirichletValues();
  rhs.Collect();
  systemMatrix.ClearDirichletValues();
}

void ILinearAssemble::AssembleSystem(Matrix &systemMatrix, Vectors &rhs) const {
  rhs = 0;
  systemMatrix = 0;
  TRY {
    for (cell ce = rhs.cells(); ce != rhs.cells_end(); ++ce)
      AssembleSystem(ce, systemMatrix, rhs);
  }
  CATCH("Error in Residual")
  rhs.ClearDirichletValues();
  rhs.Collect();
  systemMatrix.ClearDirichletValues();
}

bool isAssembleConsistent(const IAssemble &A, Vector &u) {
  // Dirichlet Boundaries are not considered in numerical derivatives
  u.ClearDirichletFlags();
  double h = 1e-6;

  Vector defect(u);
  Matrix jacobi(u);
  RMatrix numJacobi(u.size());

  A.Residual(u, defect);
  A.Jacobi(u, jacobi);

  for (int i = 0; i < u.size(); ++i) {
    Vector vp(0.0, u);
    Vector vm(0.0, u);
    vp()[i] += h;
    vm()[i] -= h;
    Vector posDefect(vp);
    Vector negDefect(vm);
    A.Residual(vp, posDefect);
    A.Residual(vm, negDefect);

    for (int j = 0; j < u.size(); ++j)
      numJacobi(j, i) = (posDefect()[j] - negDefect()[j]) / (2.0 * h);
  }

  SparseMatrix matJacobi(jacobi);

  for (int r = 0; r < matJacobi.rows(); ++r) {
    for (int c = 0; c < matJacobi.cols(); ++c) {
      if (abs(numJacobi(r, c) - matJacobi(r, c)) > h) {
        Warning("Numerical and Assembled Jacobian do not match in (row: " + std::to_string(r)
                + "|col: " + std::to_string(c)
                + ") \n numerical J(r,c) = " + std::to_string(numJacobi(r, c))
                + "\n assembled J(r,c) = " + std::to_string(matJacobi(r, c)))

                mout
            << numJacobi << endl
            << endl
            << matJacobi << endl;
        return false;
      }
    }
  }
  return true;
}
