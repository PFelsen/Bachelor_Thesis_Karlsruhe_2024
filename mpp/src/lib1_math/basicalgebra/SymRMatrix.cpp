#include "SymRMatrix.hpp"
#include <cmath>
#include "Parallel.hpp"
#include "SaveLoad.hpp"

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
// factorization of symmetric matrix
extern "C" void dsytrf_(char *UPLO, int *N, void *A, int *LDA, int *IPIV, void *WORK, void *LWORK,
                        int *INFO);
// Invert matrix using factorization
extern "C" void dsytri_(char *UPLO, int *N, void *A, int *LDA, int *IPIV, void *WORK, int *INFO);

//==================================================================================================

template<typename REAL>
SymRMatrixT<REAL>::SymRMatrixT(int dim) : dim(dim), z((dim * (dim + 1)) / 2, REAL{}) {}

template<typename REAL>
void SymRMatrixT<REAL>::resize(int dim) {
  z.resize((dim * (dim + 1)) / 2);
  this->dim = dim;
}

template<typename REAL>
void SymRMatrixT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
SymRMatrixT<REAL> &SymRMatrixT<REAL>::Identity() {
  for (int i = 0; i < dim; i++)
    for (int j = i; j < dim; j++)
      (*this)(i, j) = REAL(i == j);
  return *this;
}

template<>
SymRMatrixT<double> &SymRMatrixT<double>::Invert() {
  double A[dim * dim];
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      A[i * dim + j] = (*this)(i, j);

  int I[dim];
  int info;
  int lwork = 2 * dim;
  double work[2 * lwork];
  char uplo = 'L';
  dsytrf_(&uplo, &dim, A, &dim, I, work, &lwork, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }
  dsytri_(&uplo, &dim, A, &dim, I, work, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }

  for (int i = 0; i < dim; ++i)
    for (int j = i; j < dim; ++j)
      (*this)(i, j) = A[i * dim + j];

  return *this;
}

template<typename REAL>
Saver &SymRMatrixT<REAL>::save(Saver &saver) const {
  saver << Dim();
  for (int i = 0; i < z.size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &SymRMatrixT<REAL>::load(Loader &loader) {
  int N;
  loader >> N;
  resize(N);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template class SymRMatrixT<double>;

#ifdef BUILD_IA

template class SymRMatrixT<IAInterval>;

SymRMatrix mid(const IASymRMatrix &IA_A) {
  SymRMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A(mid(IA_A(i, j)), i, j);
  return A;
}

SymRMatrix sup(const IASymRMatrix &IA_A) {
  SymRMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A(sup(IA_A(i, j)), i, j);
  return A;
}

SymRMatrix inf(const IASymRMatrix &IA_A) {
  SymRMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A(inf(IA_A(i, j)), i, j);
  return A;
}

#endif // BUILD_IA
