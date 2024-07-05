#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include "CMatrix.hpp"
#include "RMatrix.hpp"

/// Computes the Cholesky decomposition of a real symmetric matrix A=L*L^T
template<typename T>
RMatrixT<T> Cholesky(const SymRMatrixT<T> &A) {
  int R = A.Dim();
  RMatrixT<T> L(R);
  for (int i = 0; i < R; ++i) {
    T sum = A(i, i);
    for (int k = 0; k < i; ++k)
      sum -= L(i, k) * L(i, k);
    L[i][i] = sqrt(sum);
    if (abs(L[i][i]) < GeometricTolerance) THROW("singular matrix in Cholesky")
    for (int j = i + 1; j < R; j++) {
      T sum = A(i, j);
      for (int k = 0; k < i; k++)
        sum -= L(j, k) * L(i, k);
      L[j][i] = sum / L[i][i];
    }
  }
  return L;
}

/// Computes the Cholesky decomposition of a real symmetric matrix A=L*L^T
template<typename T>
RMatrixT<T> Cholesky(const RMatrixT<T> &A) {
  if (A.rows() != A.cols()) THROW("non square matrix in Cholesky")
  SymRMatrixT<T> A_sym(A.rows());
  for (int i = 0; i < A.rows(); ++i) {
    A_sym(i, i) = A[i][i];
    for (int j = 0; j < i; ++j) {
      if (abs(A[i][j] - A[j][i]) > GeometricTolerance) THROW("non symmtetric matrix in Cholesky")
      A_sym(A[i][j], i, j);
    }
  }
  return Cholesky(A_sym);
}

/// Computes the Cholesky decomposition of a complex hermitian matrix A=L*L^H
template<typename T>
CMatrixT<T> Cholesky(const HermCMatrixT<T> &A) {
  int R = A.Dim();
  CMatrixT<T> L(R);
  for (int i = 0; i < R; ++i) {
    T sum = std::real(A(i, i));
    for (int k = 0; k < i; ++k)
      sum -= std::norm(L(i, k));
    L[i][i] = sqrt(sum);
    if (abs(L[i][i]) < GeometricTolerance) THROW("singular matrix in Cholesky")
    for (int j = i + 1; j < R; j++) {
      std::complex<T> sum = std::conj(A(i, j));
      for (int k = 0; k < i; k++)
        sum -= L(j, k) * std::conj(L(i, k));
      L[j][i] = sum / L[i][i];
    }
  }
  return L;
}

/// Computes the Cholesky decomposition of a complex hermitian matrix A=L*L^H
template<typename T>
CMatrixT<T> Cholesky(const CMatrixT<T> &A) {
  if (A.rows() != A.cols()) THROW("non square matrix in Cholesky")
  HermCMatrixT<T> A_herm(A.rows());
  for (int i = 0; i < A.rows(); ++i) {
    if (std::imag(A[i][i]) != T(0.0)) THROW("non real diagonal entry in Cholesky")
    A_herm(A[i][i], i, i);
    for (int j = 0; j < i; ++j) {
      if (abs(A[i][j] - std::conj(A[j][i])) > GeometricTolerance)
        Exit("non hermitian matrix in Cholesky") A_herm(A[i][j], i, j);
    }
  }
  return Cholesky(A_herm);
}

namespace mpp_cholesky_solve {

/// Solves the linear system L*L^T*x=b, where b needs to be stored in x and L is a lower triagular
/// matrix
template<typename T>
void forward_backward(int dim, const RMatrixT<T> &L, RVectorT<T> &x) {
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j)
      x[i] -= L[i][j] * x[j];
    x[i] /= L[i][i];
  }
  for (int i = dim - 1; i >= 0; --i) {
    for (int j = i + 1; j < dim; ++j)
      x[i] -= L[j][i] * x[j];
    x[i] /= L[i][i];
  }
}

/// Solves the linear system L*L^H*x=b, where b needs to be stored in x and L is a lower triagular
/// matrix
template<typename T>
void forward_backward(int dim, const CMatrixT<T> &L, CVectorT<T> &x) {
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j)
      x[i] -= L[i][j] * x[j];
    x[i] /= L[i][i];
  }
  for (int i = dim - 1; i >= 0; --i) {
    for (int j = i + 1; j < dim; ++j)
      x[i] -= std::conj(L[j][i]) * x[j];
    x[i] /= L[i][i];
  }
}
} // namespace mpp_cholesky_solve

/// Solves the linear system A*x=b for a real symmetric matrix A via Cholesky decomposition
template<typename T>
RVectorT<T> linearSystem(const SymRMatrixT<T> &A, const RVectorT<T> &b) {
  if (A.Dim() != b.size()) THROW("Dimensions do not fit in linearSystem")
  RMatrixT<T> L = Cholesky(A);
  RVectorT<T> x(b);
  mpp_cholesky_solve::forward_backward(A.Dim(), L, x);
  return x;
}

/// Solves the linear System A*x=b for a real symmetric matrix A via Cholesky decomposition
template<typename T>
RVectorT<T> linearSystem(const RMatrixT<T> &A, const RVectorT<T> &b) {
  if (A.cols() != b.size()) THROW("Dimensions do not fit in linearSystem")
  RMatrixT<T> L = Cholesky(A);
  RVectorT<T> x(b);
  mpp_cholesky_solve::forward_backward(A.cols(), L, x);
  return x;
}

/// Solves the linear system A*x=b for a complex hermitian matrix A via Cholesky decomposition
template<typename T>
CVectorT<T> linearSystem(const HermCMatrixT<T> &A, const CVectorT<T> &b) {
  if (A.Dim() != b.size()) THROW("Dimensions do not fit in linearSystem")
  CMatrixT<T> L = Cholesky(A);
  CVectorT<T> x(b);
  mpp_cholesky_solve::forward_backward(A.Dim(), L, x);
  return x;
}

/// Solves the linear system A*x=b for a complex hermitian matrix A via Cholesky decomposition
template<typename T>
CVectorT<T> linearSystem(const CMatrixT<T> &A, const CVectorT<T> &b) {
  if (A.cols() != b.size()) THROW("Dimensions do not fit in linearSystem")
  CMatrixT<T> L = Cholesky(A);
  CVectorT<T> x(b);
  mpp_cholesky_solve::forward_backward(A.cols(), L, x);
  return x;
}

#endif // CHOLESKY_HPP
