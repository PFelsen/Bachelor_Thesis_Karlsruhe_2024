#include "Spectrum.hpp"

#include <cmath>

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
extern "C" void dsyev_(char *JOBZ, char *UPLO, int *N, void *A, int *LDA, void *W, void *WORK,
                       int *LWORK, int *INFO, int ch1, int ch2);

// Implement this for a selected amount of eigenvectors.
/*extern "C" void dsyevx_( char* jobz, char* range, char* uplo, int* n, double* a,
                    int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                    int* m, double* w, double* z, int* ldz, double* work, int* lwork,
                    int* iwork, int* ifail, int* info ); */

extern "C" void dsygv_(int *ITYPE, char *JOBZ, char *UPLO, int *N, void *A, int *LDA, void *B,
                       int *LDB, void *W, void *WORK, int *LWORK, int *INFO, int ch1, int ch2);
extern "C" void dgeev_(char *JOBVL, char *JOBVR, int *N, void *A, int *LDA, void *WR, void *WI,
                       void *VL, int *LDVL, void *VR, int *LDVR, void *WORK, int *LWORK, int *INFO,
                       int ch1, int ch2);
extern "C" void zheev_(char *JOBZ, char *UPLO, int *N, void *A, int *LDA, void *W, void *WORK,
                       int *LWORK, void *RWORK, int *INFO, int ch1, int ch2);
extern "C" void zhegv_(int *ITYPE, char *JOBZ, char *UPLO, int *N, void *A, int *LDA, void *B,
                       int *LDB, void *W, void *WORK, int *LWORK, void *RWORK, int *INFO, int ch1,
                       int ch2);
extern "C" void zgeev_(char *JOBVL, char *JOBVR, int *N, void *A, int *LDA, void *W, void *VL,
                       int *LDVL, void *VR, int *LDVR, void *WORK, int *LWORK, void *RWORK,
                       int *INFO, int ch1, int ch2);

//==================================================================================================

void EVrealT(const SymRMatrix &A, Eigenvalues &lambda, Eigenvectors *x, int R, char jobz) {
  char uplo = 'U';
  int n = R, lwork = 3 * n, info;
  std::vector<double> w(n, 0.0);
  std::vector<double> work(2 * lwork, 0.0);
  std::vector<double> a(n * n);

  for (int i = 0; i < R; i++)
    for (int j = 0; j < R; j++)
      a[i * n + j] = A(i, j);

  dsyev_(&jobz, &uplo, &n, a.data(), &n, w.data(), work.data(), &lwork, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: The algorithm failed to converge; " + std::to_string(info)
            + " off-diagonal elements of an intermediate tridiagonal form did not converge to zero")
    }
  }

  lambda.resize(R);
  if (x) x->resize(R);
  for (int i = 0; i < R; i++) {
    lambda[i] = w[i];
    if (x)
      for (int j = 0; j < R; j++)
        (*x)[j][i] = a[i * R + j];
  }
}

void EVreal(const SymRMatrix &A, Eigenvalues &lambda, Eigenvectors &x) {
  EVrealT(A, lambda, &x, A.Dim(), 'V');
}

void EVreal(const SymRMatrix &A, Eigenvalues &lambda) {
  Eigenvectors *x = nullptr;
  EVrealT(A, lambda, x, A.Dim(), 'N');
}

void EVrealT(const RMatrix &A, CEigenvalues &lambda, CEigenvectors *x, int R, char jobvr) {
  char jobvl = 'N';
  int n = R, lwork = 4 * n, info;
  std::vector<double> wr(n, 0.0);
  std::vector<double> wi(n, 0.0);
  std::vector<double> vl(n, 0.0);
  std::vector<double> vr(n * n, 0.0);
  std::vector<double> work(4 * lwork, 0.0);
  std::vector<double> a(n * n);

  for (int i = 0; i < R; i++)
    for (int j = 0; j < R; j++)
      a[j * n + i] = A(i, j);

  dgeev_(&jobvl, &jobvr, &n, a.data(), &n, wr.data(), wi.data(),
         vl.data(), &n, vr.data(), &n, work.data(),
         &lwork, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: The QR algorithm failed to compute all the eigenvalues, and no "
            "eigenvectors have been computed; elements "
            + std::to_string(info) + "+1:N of WR and WI contain eigenvalues which have converged. ")
    }
  }

  lambda.resize(R);
  if (x) x->resize(R);
  for (int i = 0; i < R; i++) {
    lambda[i] = std::complex<double>(wr[i], wi[i]);
    if (x) {
      bool omitNext = false;
      if (wi[i] > 0.0) { // conjugated complex eigenvalue pair (part with positive imag. part)
        for (int j = 0; j < R; j++)
          (*x)[j][i] = std::complex<double>(vr[i * R + j], vr[(i + 1) * R + j]);
      } else if (wi[i] < 0.0) { // conjugated complex eigenvalue pair (with negative imag. part)
        for (int j = 0; j < R; j++)
          (*x)[j][i] = std::complex<double>(vr[(i - 1) * R + j], -vr[i * R + j]);
      } else { // real eigenvalue
        for (int j = 0; j < R; j++)
          (*x)[j][i] = vr[i * R + j];
      }
    }
  }
}

void EVreal(const RMatrix &A, CEigenvalues &lambda, CEigenvectors &x) {
  if (A.cols() != A.rows()) THROW("No square matrix")
  EVrealT(A, lambda, &x, A.rows(), 'V');
}

void EVreal(const RMatrix &A, CEigenvalues &lambda) {
  if (A.cols() != A.rows()) THROW("No square matrix")
  CEigenvectors *x = nullptr;
  EVrealT(A, lambda, x, A.rows(), 'N');
}

void EVrealT(const SymRMatrix &A, const SymRMatrix &B, Eigenvalues &lambda, Eigenvectors *x, int R,
             char jobz) {
  char uplo = 'U';
  int itype = 1, n = R, lwork = 3 * n, info;
  std::vector<double> w(n, 0.0);
  std::vector<double> work(2 * lwork, 0.0);
  std::vector<double> a(n * n);
  std::vector<double> b(n * n);

  for (int i = 0; i < R; i++)
    for (int j = 0; j < R; j++) {
      a[i * n + j] = A(i, j);
      b[i * n + j] = B(i, j);
    }

  dsygv_(&itype, &jobz, &uplo, &n, a.data(), &n, b.data(), &n, w.data(),
         work.data(), &lwork, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      if (info <= R) {
        THROW(
            "Error in LAPACK: The algorithm failed to converge; " + std::to_string(info)
            + " off-diagonal elements of an intermediate tridiagonal form did not converge to zero")
      } else {
        THROW("Error in LAPACK: The algorithm failed to converge; the leading minor of order "
              + std::to_string(info)
              + " of B is not positive definite. The factorization of B could not be completed and "
                "no eigenvalues or eigenvectors were computed")
      }
    }
  }

  lambda.resize(R);
  if (x) x->resize(R);
  for (int i = 0; i < R; i++) {
    lambda[i] = w[i];
    if (x)
      for (int j = 0; j < R; j++)
        (*x)[j][i] = a[i * R + j];
  }
}

void EVreal(const SymRMatrix &A, const SymRMatrix &B, Eigenvalues &lambda, Eigenvectors &x) {
  EVrealT(A, B, lambda, &x, A.Dim(), 'V');
}

void EVreal(const SymRMatrix &A, const SymRMatrix &B, Eigenvalues &lambda) {
  Eigenvectors *x = nullptr;
  EVrealT(A, B, lambda, x, A.Dim(), 'N');
}

void EVcomplexT(const HermCMatrix &A, Eigenvalues &lambda, CEigenvectors *x, int R, char jobz) {
  char uplo = 'U';
  int itype = 1, n = R, lwork = 2 * n, info;
  std::vector<double> w(n, 0.0);
  std::vector<double> work(2 * lwork, 0.0);
  std::vector<double> rwork(3 * n, 0.0);
  std::vector<std::complex<double>> a(n * n);

  for (int i = 0; i < R; i++)
    for (int j = 0; j < R; j++)
      a[i * n + j] = A(i, j);

  zheev_(&jobz, &uplo, &n, a.data(), &n, w.data(), work.data(), &lwork,
         rwork.data(), &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: The algorithm failed to converge; " + std::to_string(info)
            + " off-diagonal elements of an intermediate tridiagonal form did not converge to zero")
    }
  }

  lambda.resize(R);
  if (x) x->resize(R);
  for (int i = 0; i < R; i++) {
    lambda[i] = w[i];
    if (x)
      for (int j = 0; j < R; j++)
        (*x)[j][i] = a[i * R + j];
  }
}

void EVcomplex(const HermCMatrix &A, Eigenvalues &lambda, CEigenvectors &x) {
  EVcomplexT(A, lambda, &x, A.Dim(), 'V');
}

void EVcomplex(const HermCMatrix &A, Eigenvalues &lambda) {
  CEigenvectors *x = nullptr;
  EVcomplexT(A, lambda, x, A.Dim(), 'N');
}

void EVcomplexT(const CMatrix &A, CEigenvalues &lambda, CEigenvectors *x, int R, char jobvr) {
  char jobvl = 'N';
  int n = R, lwork = 4 * n, info;
  std::vector<std::complex<double>> w(n, 0.0);
  std::vector<std::complex<double>> vl(n, 0.0);
  std::vector<std::complex<double>> vr(n * n, 0.0);
  std::vector<std::complex<double>> work(4 * lwork, 0.0);
  std::vector<double> workr(2 * n, 0.0);
  std::vector<std::complex<double>> a(n * n);

  for (int i = 0; i < R; i++)
    for (int j = 0; j < R; j++)
      a[j * n + i] = A(i, j);


  zgeev_(&jobvl, &jobvr, &n, a.data(), &n, w.data(), vl.data(),
         &n, vr.data(), &n, work.data(), &lwork,
         workr.data(), &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      THROW("Error in LAPACK: The QR algorithm failed to compute all the eigenvalues, and no "
            "eigenvectors have been computed; elements "
            + std::to_string(info) + "+1:N of WR and WI contain eigenvalues which have converged. ")
    }
  }

  lambda.resize(R);
  if (x) x->resize(R);
  for (int i = 0; i < R; i++) {
    lambda[i] = w[i];
    if (x)
      for (int j = 0; j < R; j++)
        (*x)[j][i] = vr[i * R + j];
  }
}

void EVcomplex(const CMatrix &A, CEigenvalues &lambda, CEigenvectors &x) {
  if (A.cols() != A.rows()) THROW("No square matrix")
  EVcomplexT(A, lambda, &x, A.rows(), 'V');
}

void EVcomplex(const CMatrix &A, CEigenvalues &lambda) {
  if (A.cols() != A.rows()) THROW("No square matrix")
  CEigenvectors *x = nullptr;
  EVcomplexT(A, lambda, x, A.rows(), 'N');
}

void EVcomplexT(const HermCMatrix &A, const HermCMatrix &B, Eigenvalues &lambda, CEigenvectors *x,
                int R, char jobz) {
  char uplo = 'U';
  int itype = 1, n = R, lwork = 2 * n, info;
  std::vector<double> w(n, 0.0);
  std::vector<double> work(2 * lwork, 0.0);
  std::vector<double> rwork(3 * n, 0.0);

  std::vector<std::complex<double>> a(n * n);
  std::vector<std::complex<double>> b(n * n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      a[i * n + j] = A(i, j);
      b[i * n + j] = B(i, j);
    }

  zhegv_(&itype, &jobz, &uplo, &n, a.data(), &n, b.data(), &n, w.data(),
         work.data(), &lwork, rwork.data(), &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    } else if (info > 0) {
      if (info <= R) {
        THROW(
            "Error in LAPACK: The algorithm failed to converge; " + std::to_string(info)
            + " off-diagonal elements of an intermediate tridiagonal form did not converge to zero")
      } else {
        THROW("Error in LAPACK: The algorithm failed to converge; the leading minor of order "
              + std::to_string(info)
              + " of B is not positive definite. The factorization of B could not be completed and "
                "no eigenvalues or eigenvectors were computed")
      }
    }
  }

  lambda.resize(R);
  if (x) x->resize(R);
  for (int i = 0; i < R; i++) {
    lambda[i] = w[i];
    if (x)
      for (int j = 0; j < R; j++)
        (*x)[j][i] = a[i * R + j];
  }
}

void EVcomplex(const HermCMatrix &A, const HermCMatrix &B, Eigenvalues &lambda, CEigenvectors &x) {
  EVcomplexT(A, B, lambda, &x, A.Dim(), 'V');
}

void EVcomplex(const HermCMatrix &A, const HermCMatrix &B, Eigenvalues &lambda) {
  CEigenvectors *x = nullptr;
  EVcomplexT(A, B, lambda, x, A.Dim(), 'N');
}

Spectrum::Spectrum(const SymRMatrix &A, bool computeEVec) :
    computeEVec(computeEVec), lambda(A.Dim()) {
  if (computeEVec) {
    e1 = std::make_unique<Eigenvectors>(A.Dim());
    EVreal(A, lambda, *e1);
  } else {
    EVreal(A, lambda);
  }
}

Spectrum::Spectrum(const HermCMatrix &A, bool computeEVec) :
    computeEVec(computeEVec), lambda(A.Dim()) {
  if (computeEVec) {
    e2 = std::make_unique<CEigenvectors>(A.Dim());
    EVcomplex(A, lambda, *e2);
  } else {
    EVcomplex(A, lambda);
  }
}

Spectrum::Spectrum(const SymRMatrix &A, const SymRMatrix &B, bool computeEVec) :
    computeEVec(computeEVec), lambda(A.Dim()) {
  if (computeEVec) {
    e1 = std::make_unique<Eigenvectors>(A.Dim());
    EVreal(A, B, lambda, *e1);
  } else {
    EVreal(A, B, lambda);
  }
}

Spectrum::Spectrum(const HermCMatrix &A, const HermCMatrix &B, bool computeEVec) :
    computeEVec(computeEVec), lambda(A.Dim()) {
  if (computeEVec) {
    e2 = std::make_unique<CEigenvectors>(A.Dim());
    EVcomplex(A, B, lambda, *e2);
  } else {
    EVcomplex(A, B, lambda);
  }
}

Eigenvectors &Spectrum::getEigenvectors() const {
  if (!e1) THROW("No such eigenvector")
  return *e1;
}

CEigenvectors &Spectrum::getCEigenvectors() const {
  if (!e2) THROW("No such eigenvector")
  return *e2;
}

Eigenvalue Spectrum::operator()(int i, Eigenvector &w) const {
  if (!e1) THROW("No such eigenvector")
  w.resize(lambda.size());
  for (int j = 0; j < lambda.size(); ++j)
    w[j] = (*e1)[j][i];
  return lambda[i];
}

Eigenvalue Spectrum::operator()(int i, CEigenvector &w) const {
  if (!e1) THROW("No such eigenvector")
  w.resize(lambda.size());
  for (int j = 0; j < lambda.size(); ++j)
    w[j] = (*e2)[j][i];
  return lambda[i];
}

Eigenvalue Spectrum::min() const {
  int i = 0;
  while (std::abs(lambda[i] - 1.0) < 1e-8)
    ++i;
  return lambda[i];
}

double Spectrum::absmin() const {
  double s = infty;
  for (int j = 0; j < lambda.size(); ++j)
    if (std::abs(lambda[j] - 1.0) > 1e-8)
      if (s > std::abs(lambda[j])) s = std::abs(lambda[j]);
  return s;
}

Eigenvalue Spectrum::max() const { return lambda[lambda.size() - 1]; }

double Spectrum::absmax() const {
  double s = 0;
  for (int j = 0; j < lambda.size(); ++j)
    if (std::abs(lambda[j] - 1.0) > 1e-8) s = std::max(std::abs(lambda[j]), s);
  return s;
}

double Spectrum::cond() const { return absmax() / absmin(); }