#include "QRDecomposition.hpp"

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
extern "C" void dgeqrf_(int *M, int *N, void *A, int *LDA, void *TAU, void *WORK, int *LWORK,
                        int *INFO, int ch1, int ch2);

//==================================================================================================

template<>
QRMatricesT<double> QRDecomposition(const RMatrix &A) {
  QRMatricesT<double> QR(A);
  int lwork = 3 * QR.N, info;
  double work[2 * lwork];

  dgeqrf_(&QR.M, &QR.N, QR.A.data(), &QR.M, QR.TAU.data(), work, &lwork, &info, 1, 1);

  if (info != 0) {
    if (info < 0) {
      THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
    }
  }
  return QR;
}