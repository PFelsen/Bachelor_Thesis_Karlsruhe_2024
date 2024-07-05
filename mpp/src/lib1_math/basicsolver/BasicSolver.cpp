#include "BasicSolver.hpp"

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
extern "C" void dgesv_(int *N, int *Nrhs, void *A, int *LDA, int *IPIV, void *B, int *LDB,
                       int *INFO, int ch1, int ch2);
extern "C" void dsysv_(char *UPLO, int *N, int *Nrhs, void *A, int *LDA, int *IPIV, void *B,
                       int *LDB, void *WORK, int *LWORK, int *INFO, int ch1, int ch2);
extern "C" void zgesv_(int *N, int *Nrhs, void *A, int *LDA, int *IPIV, void *B, int *LDB,
                       int *INFO, int ch1, int ch2);
extern "C" void zhesv_(char *UPLO, int *N, int *Nrhs, void *A, int *LDA, int *IPIV, void *B,
                       int *LDB, void *WORK, int *LWORK, int *INFO, int ch1, int ch2);

//==================================================================================================

void LAPACKSolve(const RMatrix &A, std::vector<double> &b, int Nb) {
  if (A.rows() != A.cols()) THROW("Not a square matrix!")

  int n = A.rows(), info;
  double a[n * n];
  int piv[n];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i + j * n] = A(i, j);

  dgesv_(&n, &Nb, a, &n, piv, b.data(), &n, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: U(" + std::to_string(info) + ", " + std::to_string(info) + ")"
            + " is exactly zero. The factorization has been completed, but the factor U is exactly "
            + " singular, so the solution could not be computed.")
    }
  }
}

void LAPACKSolve(const SymRMatrix &A, std::vector<double> &b, int Nb) {
  char uplo = 'U';
  int n = A.rows(), lwork = 3 * n, info;
  double work[2 * lwork];
  double a[n * n];
  int piv[n];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i * n + j] = A(i, j);

  dsysv_(&uplo, &n, &Nb, a, &n, piv, b.data(), &n, work, &lwork, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: D(" + std::to_string(info) + ", " + std::to_string(info) + ")"
            + " is exactly zero. The factorization has been completed, but the block diagonal"
            + " matrix D is exactly singular, so the solution could not be computed.")
    }
  }
}

void LAPACKSolve(const AntisymRMatrix &A, std::vector<double> &b, int Nb) {
  int n = A.rows(), info;
  double a[n * n];
  int piv[n];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i + j * n] = A(i, j);

  dgesv_(&n, &Nb, a, &n, piv, b.data(), &n, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: U(" + std::to_string(info) + ", " + std::to_string(info) + ")"
            + " is exactly zero. The factorization has been completed, but the factor U is exactly "
            + " singular, so the solution could not be computed.")
    }
  }
}

void LAPACKSolve(const CMatrix &A, std::vector<std::complex<double>> &b, int Nb) {
  if (A.rows() != A.cols()) THROW("Not a square matrix!")

  int n = A.rows(), info;
  std::complex<double> a[n * n];
  int piv[n];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i + j * n] = A(i, j);

  zgesv_(&n, &Nb, a, &n, piv, b.data(), &n, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: U(" + std::to_string(info) + ", " + std::to_string(info) + ")"
            + " is exactly zero. The factorization has been completed, but the factor U is exactly"
            + " singular, so the solution could not be computed.")
    }
  }
}

void LAPACKSolve(const HermCMatrix &A, std::vector<std::complex<double>> &b, int Nb) {
  char uplo = 'U';
  int n = A.rows(), lwork = 3 * n, info;
  std::complex<double> work[2 * lwork];
  std::complex<double> a[n * n];
  int piv[n];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i * n + j] = A(i, j);

  zhesv_(&uplo, &n, &Nb, a, &n, piv, b.data(), &n, work, &lwork, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: D(" + std::to_string(info) + ", " + std::to_string(info) + ")"
            + " is exactly zero. The factorization has been completed, but the block diagonal"
            + " matrix D is exactly singular, so the solution could not be computed.")
    }
  }
}