#include "LOBPCG.hpp"

void LOBPCG::ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                            Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) {
  Vectors &a1 = additional[0];

  setMatricesSmall(u, a1, A, B, tmp1, tmp2);
  EVreal(*a, *b, lambda, *e);

  Vectors u_copy(u);
  for (int r = 0; r < u.size(); r++) {
    u[r] = 0;
    for (int s = 0, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted) {
      u[r] += (*e)[s][r] * u_copy[s];
      u[r] += (*e)[s_shifted][r] * a1[s];
    }
  }
}

void LOBPCG::update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A,
                    Matrix &B, Vector &residual,
                    std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
                    double defect, int r) {
  Vectors &a1 = additional[0];

  if (defect < eps) {
    randVector(a1[r], B, *solver, residual);
    ++evConverged;
  } else {
    a1[r] = (*solver) * residual;
  }
}

void LOBPCGExtended::ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional,
                                    Eigenvalues &lambda, Operator &A, Operator &B, Vector &tmp1,
                                    Vector &tmp2) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  if (firstStep) {
    setMatricesSmall(u, a1, A, B, tmp1, tmp2);
  } else {
    setMatricesLarge(u, a1, a2, A, B, tmp1, tmp2);
  }

  EVreal(*a, *b, lambda, *e);

  Eigenfcts u_copy(u);
  if (firstStep) {
    for (int r = 0; r < u.size(); r++) {
      a2[r] = 0;
      for (int s = 0, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted)
        a2[r] += (*e)[s_shifted][r] * a1[s];
      u[r] = a2[r];
      for (int s = 0; s < u.size(); ++s)
        u[r] += (*e)[s][r] * u_copy[s];
    }
    a = std::make_unique<SymRMatrix>(3 * u.size());
    b = std::make_unique<SymRMatrix>(3 * u.size());
    e = std::make_unique<RMatrix>(3 * u.size());
    firstStep = false;
  } else {
    Vectors a2_copy(a2);
    for (int r = 0; r < u.size(); r++) {
      a2[r] = 0;
      for (int s = 0, s_shifted = u.size(), s_shifted2 = 2 * s_shifted; s < u.size();
           ++s, ++s_shifted, ++s_shifted2) {
        a2[r] += (*e)[s_shifted][r] * a1[s];
        a2[r] += (*e)[s_shifted2][r] * a2_copy[s];
      }
      u[r] = a2[r];
      for (int s = 0; s < u.size(); ++s)
        u[r] += (*e)[s][r] * u_copy[s];
    }
  }
}

void LOBPCGExtended::update(
    Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
    Vector &residual, std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
    double defect, int r) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  if (defect < eps) {
    randVector(a1[r], B, *solver, residual);
    randVector(a2[r], B, *solver, residual);
    ++evConverged;
  } else {
    a1[r] = (*solver) * residual;
  }
}
