#include "LOBPCGSelective.hpp"

void LOBPCGSelective::ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional,
                                     Eigenvalues &lambda, Operator &A, Operator &B, Vector &tmp1,
                                     Vector &tmp2) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  if (firstStep) {
    setMatricesSmallSelective(u, a1, A, B, tmp1, tmp2);
  } else {
    setMatricesLargeSelective(u, a1, a2, A, B, tmp1, tmp2);
  }

  EVreal(*a, *b, lambda, *e);

  Eigenfcts u_copy(u);
  if (firstStep) {
    for (int r = 0; r < u.size(); r++) {
      if (selection[r] == 2) {
        a2[r] = 0;
        for (int s = 0, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted)
          a2[r] += (*e)[s_shifted][r] * a1[s];
        u[r] = a2[r];
        for (int s = 0; s < u.size(); ++s)
          u[r] += (*e)[s][r] * u_copy[s];
      }
    }
    a->resize(3 * u.size());
    b->resize(3 * u.size());
    e->resize(3 * u.size());
    for (int r = 0; r < u.size(); ++r)
      selection[r] = (selection[r] == 1) ? 3 : 2;
    firstStep = false;
  } else {
    Vectors a2_copy(a2);
    for (int r = 0; r < u.size(); r++) {
      if (selection[r] == 2) {
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
}

bool LOBPCGSelective::update(
    Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
    Vector &residual, Vector &tmp,
    std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    int fev) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  evConverged = 0;
  for (int r = 0; r < u.size(); ++r) {
    if (selection[r] == 2) {
      double defect = ritzDefect(u[r], residual, lambda[r], A, B, tmp);
      if (defect < eps) {
        randVector(a1[r], B, *solver, residual);
        randVector(a2[r], B, *solver, residual);
        selection[r] = 1;
      }
    } else {
      selection[r] = (selection[r] == 3) ? 1 : 0;
    }
    if (selection[r] == 2) {
      a1[r] = (*solver) * residual;
    } else if (r < fev) {
      ++evConverged;
      if (evConverged >= fev) return true;
    }
  }
  return false;
}

void LOBPCGSelectiveFast::ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional,
                                         Eigenvalues &lambda, Operator &A, Operator &B,
                                         Vector &tmp1, Vector &tmp2) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  if (firstStep) {
    setMatricesSmall(u, a1, A, B, tmp1, tmp2);
  } else {
    setMatricesLarge(u, a1, a2, A, B, tmp1, tmp2);
  }

  EVreal(*a, *b, lambda, *e);

  Eigenfcts u_copy(u);
  Vectors a2_copy(a2);
  if (firstStep) {
    for (int r = 0; r < u.size(); r++) {
      if (selection[r] == 2) {
        a2[r] = 0;
        for (int s = 0, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted)
          a2[r] += (*e)[s_shifted][r] * a1[s];
        u[r] = a2[r];
        for (int s = 0; s < u.size(); ++s)
          u[r] += (*e)[s][r] * u_copy[s];
      } else {
        u[r] = 0;
        for (int s = 0, s_shifted = u.size(), s_shifted2 = 2 * s_shifted; s < u.size();
             ++s, ++s_shifted, ++s_shifted2) {
          u[r] += (*e)[s][r] * u_copy[s];
          u[r] += (*e)[s_shifted][r] * a1[s];
          u[r] += (*e)[s_shifted2][r] * a2_copy[s];
        }
      }
    }
    a->resize(3 * u.size());
    b->resize(3 * u.size());
    e->resize(3 * u.size());
    for (int r = 0; r < u.size(); ++r)
      selection[r] = (selection[r] == 1) ? 3 : 2;
    firstStep = false;
  } else {
    for (int r = 0; r < u.size(); r++) {
      if (selection[r] == 2) {
        a2[r] = 0;
        for (int s = 0, s_shifted = u.size(), s_shifted2 = 2 * s_shifted; s < u.size();
             ++s, ++s_shifted, ++s_shifted2) {
          a2[r] += (*e)[s_shifted][r] * a1[s];
          a2[r] += (*e)[s_shifted2][r] * a2_copy[s];
        }
        u[r] = a2[r];
        for (int s = 0; s < u.size(); ++s)
          u[r] += (*e)[s][r] * u_copy[s];
      } else {
        u[r] = 0;
        for (int s = 0, s_shifted = u.size(), s_shifted2 = 2 * s_shifted; s < u.size();
             ++s, ++s_shifted, ++s_shifted2) {
          u[r] += (*e)[s][r] * u_copy[s];
          u[r] += (*e)[s_shifted][r] * a1[s];
          u[r] += (*e)[s_shifted2][r] * a2_copy[s];
        }
      }
    }
  }
}

bool LOBPCGSelectiveFast::update(
    Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
    Vector &residual, Vector &tmp,
    std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    int fev) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  for (int r = 0; r < u.size(); ++r) {
    if (selection[r] == 2) {
      double defect = ritzDefect(u[r], residual, lambda[r], A, B, tmp);
      if (defect < eps) {
        selection[r] = 0;
        if (r < fev) {
          ++evConverged;
          if (evConverged >= fev) return true;
        }
        randVector(a1[r], B, *solver, residual);
        randVector(a2[r], B, *solver, residual);
      } else {
        a1[r] = (*solver) * residual;
      }
    }
  }
  return false;
}

void LOBPCGSelectiveVeryFast::setMatricesSmallVeryFast(Eigenfcts &u, Vectors &a1, Operator &A,
                                                       Operator &B, Vector &tmp1, Vector &tmp2) {
  for (int s = 0; s < u.size(); ++s) {
    tmp1 = A * u[s];
    tmp2 = B * u[s];
    for (int r = 0; r <= s; ++r) {
      (*a)(u[r] * tmp1, s, r);
      (*b)(u[r] * tmp2, s, r);
    }
  }
  for (int s = evConverged, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted) {
    tmp1 = A * a1[s];
    tmp2 = B * a1[s];
    for (int r = 0; r < u.size(); ++r) {
      (*a)(u[r] * tmp1, s_shifted, r);
      (*b)(u[r] * tmp2, s_shifted, r);
    }
    for (int r = evConverged, r_shifted = u.size(); r < u.size(); ++r, ++r_shifted) {
      (*a)(a1[r] * tmp1, s_shifted, r_shifted);
      (*b)(a1[r] * tmp2, s_shifted, r_shifted);
    }
  }
}

void LOBPCGSelectiveVeryFast::setMatricesLargeVeryFast(Eigenfcts &u, Vectors &a1, Vectors &a2,
                                                       Operator &A, Operator &B, Vector &tmp1,
                                                       Vector &tmp2) {
  setMatricesSmallVeryFast(u, a1, A, B, tmp1, tmp2);
  for (int s = evConverged, s_shifted = 2 * u.size() - evConverged; s < u.size();
       ++s, ++s_shifted) {
    tmp1 = A * a2[s];
    tmp2 = B * a2[s];
    for (int r = 0; r < u.size(); ++r) {
      (*a)(u[r] * tmp1, s_shifted, r);
      (*b)(u[r] * tmp2, s_shifted, r);
    }
    for (int r = evConverged, r_shifted = u.size(), r_shifted2 = 2 * u.size() - evConverged;
         r < u.size(); ++r, ++r_shifted, ++r_shifted2) {
      (*a)(a1[r] * tmp1, s_shifted, r_shifted);
      (*b)(a1[r] * tmp2, s_shifted, r_shifted);
      (*a)(a2[r] * tmp1, s_shifted, r_shifted2);
      (*b)(a2[r] * tmp2, s_shifted, r_shifted2);
    }
  }
}

void LOBPCGSelectiveVeryFast::ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional,
                                             Eigenvalues &lambda, Operator &A, Operator &B,
                                             Vector &tmp1, Vector &tmp2) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  int R = u.size();

  if (firstStep) {
    if (evConverged > 0) {
      a->resize(2 * R - evConverged);
      b->resize(2 * R - evConverged);
      e->resize(2 * R - evConverged);
    }
    setMatricesSmallVeryFast(u, a1, A, B, tmp1, tmp2);
  } else {
    int matrixSize = 3 * R - 2 * evConverged;
    if (a->Dim() != matrixSize) {
      a->resize(matrixSize);
      b->resize(matrixSize);
      e->resize(matrixSize);
    }
    setMatricesLargeVeryFast(u, a1, a2, A, B, tmp1, tmp2);
  }

  try {
    EVreal(*a, *b, lambda, *e);
  } catch (...) {
    mout << "A = " << endl << *a << endl << endl << "B = " << endl << *b << endl << endl;
    Exit("Error in ritzStepLOBPCG")
  }

  Vectors u_copy(u);
  Vectors a2_copy(a2);
  if (firstStep) {
    int r;
    for (r = 0; r < evConverged; ++r) {
      u[r] = 0;
      for (int s = 0; s < R; ++s)
        u[r] += (*e)[s][r] * u_copy[s];
      for (int s = evConverged, s_shifted = R; s < R; ++s, ++s_shifted)
        u[r] += (*e)[s_shifted][r] * a1[s];
    }
    for (; r < R; ++r) {
      a2[r] = 0;
      for (int s = evConverged, s_shifted = R; s < R; ++s, ++s_shifted)
        a2[r] += (*e)[s_shifted][r] * a1[s];
      u[r] = a2[r];
      for (int s = 0; s < R; ++s)
        u[r] += (*e)[s][r] * u_copy[s];
    }
    firstStep = false;
  } else {
    int r;
    for (r = 0; r < evConverged; ++r) {
      u[r] = 0;
      for (int s = 0; s < R; ++s) {
        u[r] += (*e)[s][r] * u_copy[s];
      }
      for (int s = evConverged, s_shifted = R, s_shifted2 = 2 * R - evConverged; s < R;
           ++s, ++s_shifted, ++s_shifted2) {
        u[r] += (*e)[s_shifted][r] * a1[s];
        u[r] += (*e)[s_shifted2][r] * a2_copy[s];
      }
    }
    for (; r < R; ++r) {
      a2[r] = 0;
      for (int s = evConverged, s_shifted = R, s_shifted2 = 2 * R - evConverged; s < R;
           ++s, ++s_shifted, ++s_shifted2) {
        a2[r] += (*e)[s_shifted][r] * a1[s];
        a2[r] += (*e)[s_shifted2][r] * a2_copy[s];
      }
      u[r] = a2[r];
      for (int s = 0; s < R; ++s)
        u[r] += (*e)[s][r] * u_copy[s];
    }
  }
}

bool LOBPCGSelectiveVeryFast::update(
    Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
    Vector &residual, Vector &tmp,
    std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    int fev) {
  Vectors &a1 = additional[0];
  Vectors &a2 = additional[1];

  for (int r = evConverged; r < u.size(); ++r) {
    double defect = ritzDefect(u[r], residual, lambda[r], A, B, tmp);
    if ((defect < eps) && (r == evConverged)) {
      evConverged++;
      if (evConverged == fev) return true;
    } else {
      a1[r] = (*solver) * residual;
    }
  }
  return false;
}
