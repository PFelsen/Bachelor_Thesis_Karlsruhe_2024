#ifndef QRDECOMPOSITION_HPP
#define QRDECOMPOSITION_HPP

#include "CMatrix.hpp"


template<typename REAL>
class QRMatricesT;

template<typename REAL>
QRMatricesT<REAL> QRDecomposition(const RMatrixT<REAL> &B);

/**
 * LAPACK QR matrices (for further details see LAPACK function dgeqrf)
 *
 * The elements on and above the diagonal of the array contain the min(M,N)-by-N upper
 * trapezoidal matrix R (R is upper triangular if m >= n); the elements below the diagonal,
 * with the array TAU, represent the orthogonal matrix Q as a product of min(m,n) elementary
 * reflectors:
 * The matrix Q is represented as a product of elementary reflectors
 *      Q = H(1) H(2) . . . H(k), where k = min(m,n).
 * Each H(i) has the form
 *      H(i) = I - tau * v * v**T
 * where tau is a real scalar, and v is a real vector with v(1:i-1) = 0 and v(i) = 1;
 * v(i+1:m) is stored on exit in A(i+1:m,i), and tau in TAU(i).
 */
template<typename REAL>
class QRMatricesT {
protected:
  int M;
  int N;
  int K;
  std::vector<REAL> A;
  std::vector<REAL> TAU;

  /// Scalar product v_j^T * y (0 <= j < M)
  REAL ProductV(int j, const RVectorT<REAL> &y) const {
    REAL sum = y[j];
    for (int i = j + 1; i < M; ++i)
      sum += A[i + j * M] * y[i];
    return sum;
  }

  REAL ProductVCol(int j, int col, const RMatrixT<REAL> &B) const {
    REAL sum = B[j][col];
    for (int i = j + 1; i < M; ++i)
      sum += A[i + j * M] * B[i][col];
    return sum;
  }

  REAL ProductVRow(int j, int row, const RMatrixT<REAL> &B) const {
    REAL sum = B[row][j];
    for (int i = j + 1; i < M; ++i)
      sum += A[i + j * M] * B[row][i];
    return sum;
  }

  /// Product H(i) * y = y - tau * v_i * (v_i^T * y)
  RVectorT<REAL> ApplyH(int j, const RVectorT<REAL> &y) const {
    RVectorT<REAL> z(y);
    REAL p = ProductV(j, y) * TAU[j];
    z[j] -= p;
    for (int i = j + 1; i < M; ++i)
      z[i] -= p * A[i + j * M];
    return z;
  }

  /// Product H(i) * B = B - tau * v_i * (v_i^T * B)
  RMatrixT<REAL> ApplyH(int j, const RMatrixT<REAL> &B) const {
    RMatrixT<REAL> Z(B);
    for (int k = 0; k < B.cols(); ++k) {
      REAL p = ProductVCol(j, k, B) * TAU[j];
      Z[j][k] -= p;
      for (int i = j + 1; i < M; ++i)
        Z[i][k] -= p * A[i + j * M];
    }
    return Z;
  }

  /// Product B * H(i) = B - tau * (B * v_i) * v_i^T
  RMatrixT<REAL> ApplyHtransposed(int j, const RMatrixT<REAL> &B) const {
    RMatrixT<REAL> Z(B);
    for (int k = 0; k < B.rows(); ++k) {
      REAL p = ProductVRow(j, k, B) * TAU[j];
      Z[k][j] -= p;
      for (int i = j + 1; i < M; ++i)
        Z[k][i] -= p * A[i + j * M];
    }
    return Z;
  }

  QRMatricesT(int n, int m) : N(n), M(m), K(std::min(N, M)), A(M * N), TAU(K) {}

  QRMatricesT(const RMatrixT<REAL> &B) :
      M(B.rows()), N(B.cols()), K(std::min(N, M)), A(M * N), TAU(K) {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        A[i + j * M] = B(i, j);
  }
public:
  /// Product Q * y
  RVectorT<REAL> ApplyQ(const RVectorT<REAL> &y) const {
    if (M != y.size()) THROW("Size does not fit")
    RVectorT<REAL> z(y);
    for (int i = K - 1; i >= 0; --i)
      z = ApplyH(i, z);
    return z;
  }

  /// Product Q * y
  RMatrixT<REAL> ApplyQ(const RMatrixT<REAL> &B) const {
    if (M != B.rows()) THROW("Size does not fit")
    RMatrixT<REAL> Z(B);
    for (int i = K - 1; i >= 0; --i)
      Z = ApplyH(i, Z);
    return Z;
  }

  /// Product y^T * Q = (Q^T * y)^T
  RVectorT<REAL> ApplyQtransposed(const RVectorT<REAL> &y) const {
    if (M != y.size()) THROW("Size does not fit")
    RVectorT<REAL> z(y);
    for (int i = 0; i < K; ++i)
      z = ApplyH(i, z);
    return z;
  }

  /// Product y^T * Q = (Q^T * y)^T
  RMatrixT<REAL> ApplyQtransposed(const RMatrixT<REAL> &B) const {
    if (M != B.cols()) THROW("Size does not fit")
    RVectorT<REAL> Z(B);
    for (int i = 0; i < K; ++i)
      Z = ApplyH(i, Z);
    return Z;
  }

  /// Product R * y
  RVectorT<REAL> ApplyR(const RVectorT<REAL> &y) const {
    if (N != y.size()) THROW("Size does not fit")
    RVectorT<REAL> z(M);
    for (int i = 0; i < M; ++i)
      for (int j = i; j < N; ++j)
        z[i] += A[i + j * M] * y[j];
    return z;
  }

  /// Product R * B
  RMatrixT<REAL> ApplyR(const RMatrixT<REAL> &B) const {
    if (N != B.rows()) THROW("Size does not fit")
    RMatrixT<REAL> Z(M, B.cols());
    for (int i = 0; i < M; ++i)
      for (int k = 0; k < B.cols(); ++k)
        for (int j = i; j < N; ++j)
          Z[i][k] += A[i + j * M] * B[j][k];
    return Z;
  }

  /// Product y^T * R = (R^T * y)^T
  RVectorT<REAL> ApplyRtransposed(const RVectorT<REAL> &y) const {
    if (M != y.size()) THROW("Size does not fit")
    RVectorT<REAL> z(N);
    for (int i = 0; i < N; ++i)
      for (int j = 0; j <= i; ++j)
        z[i] += A[j + i * M] * y[j];
    return z;
  }

  RMatrixT<REAL> R() const {
    RMatrixT<REAL> r(M, N);
    for (int i = 0; i < M; i++)
      for (int j = i; j < N; j++)
        r[i][j] = A[i + j * M];
    return r;
  }

  friend QRMatricesT QRDecomposition<REAL>(const RMatrixT<REAL> &B);
};

template<typename REAL>
inline RVectorT<REAL> operator*(const QRMatricesT<REAL> &QR, const RVectorT<REAL> &y) {
  return QR.ApplyQ(y);
}

template<typename REAL>
inline RMatrixT<REAL> operator*(const QRMatricesT<REAL> &QR, const RMatrixT<REAL> &B) {
  return QR.ApplyQ(B);
}

template<typename REAL>
inline RVectorT<REAL> operator*(const RVectorT<REAL> &y, const QRMatricesT<REAL> &QR) {
  return QR.ApplyQtransposed(y);
}

template<typename REAL>
inline RMatrixT<REAL> operator*(const RMatrixT<REAL> &B, const QRMatricesT<REAL> &QR) {
  return QR.ApplyQtransposed(B);
}

#endif // QRDECOMPOSITION_HPP
