#include "SSOR.hpp"

#include "Matrix.hpp"
#include "Vector.hpp"

void SSOR::Construct(const Matrix &_A) {
  AA = &_A;
  Matrix *A = new Matrix(_A);
  A->Accumulate();
  N = int(_A.nR()); // CONVERSION!
  for (int i = 0; i < N; ++i) {
    int d = _A.Diag(i);
    int n = _A.Dof(i);
    D.push_back(RMatrix(n));
    for (int k = 0; k < n; ++k) {
      for (int l = 0; l < n; ++l) {
        D[i][k][l] = (*A)(d)[k * n + l];
      }
      // regularization of the diagonal blocks
      if (D[i][k][k] == Scalar(0)) D[i][k][k] = 1.0;
    }
  }
  delete A;
  for (int i = 0; i < N; ++i)
    D[i].Invert();
}

void SSOR::multiply(Vector &u, const Vector &b) const {
  u.SetAccumulateFlag(false);
  // comments refer to notation of Hanke Bourgeois p.84
  Vector _b(b);
  for (int i = 0; i < N; ++i) {
    int d = AA->Diag(i);
    int n = AA->Dof(i);
    for (++d; d < AA->Diag(i + 1); ++d) {
      int j = AA->Column(d);
      int m = AA->Dof(j);
      for (int k = 0; k < n; ++k) {
        for (int l = 0; l < m; ++l) {
          // this is a_ij x^(k+1/2) but a_ij is a matrix not a scalar
          _b(i, k) -= (*AA)(d)[k * m + l] * u(j, l);
        }
      }
    }
    for (int k = 0; k < n; ++k) {
      u(i, k) = 0.0;
      for (int l = 0; l < n; ++l) {
        // this updates x^(k+1/2)
        u(i, k) += D[i][k][l] * _b(i, l);
      }
    }
    for (int k = 0; k < n; ++k)
      u(i, k) *= omega;
  }
  _b.Collect();
  for (int i = N - 1; i >= 0; --i) {
    int d = AA->Diag(i);
    int n = AA->Dof(i);
    for (int k = 0; k < n; ++k) {
      u(i, k) = 0.0;
      for (int l = 0; l < n; ++l)
        u(i, k) += D[i][k][l] * _b(i, l);
    }
    for (int k = 0; k < n; ++k)
      u(i, k) *= omega;
    for (++d; d < AA->Diag(i + 1); ++d) {
      int j = AA->Column(d);
      int m = AA->Dof(j);
      for (int k = 0; k < m; ++k) {
        for (int l = 0; l < n; ++l) {
          _b(j, k) -= (*AA)(d)[n * m + k * n + l] * u(i, l);
        }
      }
    }
  }
  u.Accumulate();
}
