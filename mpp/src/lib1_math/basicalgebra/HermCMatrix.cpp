#include "HermCMatrix.hpp"
#include "CFunctions.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
// factorization of symmetric matrix
extern "C" void zsytrf_(char *UPLO, int *N, void *A, int *LDA, int *IPIV, void *WORK, void *LWORK,
                        int *INFO);
// Invert matrix using factorization
extern "C" void zsytri_(char *UPLO, int *N, void *A, int *LDA, int *IPIV, void *WORK, int *INFO);

//==================================================================================================

template<typename REAL>
HermCMatrixT<REAL>::HermCMatrixT(int dim) : dim(dim), z((dim * (dim + 1)) / 2, COMPLEX{}) {}

template<typename REAL>
COMPLEX_TYPE<REAL> HermCMatrixT<REAL>::operator()(int i, int j) {
  if (i >= j) return z[(i * (i + 1)) / 2 + j];
  else return baConj(z[(j * (j + 1)) / 2 + i]);
}

template<typename REAL>
COMPLEX_TYPE<REAL> HermCMatrixT<REAL>::operator()(int i, int j) const {
  if (i >= j) return z[(i * (i + 1)) / 2 + j];
  else return baConj(z[(j * (j + 1)) / 2 + i]);
}

template<typename REAL>
void HermCMatrixT<REAL>::operator()(const COMPLEX &a, int i, int j) {
  if (i == j) {
    if (baImag(a) == REAL{}) z[(i * (i + 1)) / 2 + j] = a;
    else THROW("Diagonal must be real!")
  } else if (i > j) z[(i * (i + 1)) / 2 + j] = a;
  else if (i < j) z[(j * (j + 1)) / 2 + i] = baConj(a);
}

template<typename REAL>
void HermCMatrixT<REAL>::resize(int dim) {
  z.resize((dim * (dim + 1)) / 2);
  this->dim = dim;
}

template<typename REAL>
HermCMatrixT<REAL> &HermCMatrixT<REAL>::conj() {
  for (int i = 0; i < z.size(); ++i)
    z[i] = baConj(z[i]);
  return *this;
}

template<typename REAL>
SymRMatrixT<REAL> HermCMatrixT<REAL>::real() const {
  SymRMatrixT<REAL> re(dim);
  for (int i = 0; i < z.size(); ++i)
    re[i] = baReal(z[i]);
  return re;
}

template<typename REAL>
AntisymRMatrixT<REAL> HermCMatrixT<REAL>::imag() const {
  AntisymRMatrixT<REAL> im(dim);
  for (int i = 0; i < dim; ++i)
    for (int j = i + 1; j < dim; ++j)
      im(baImag((*this)(i, j)), i, j);
  return im;
}

template<typename REAL>
void HermCMatrixT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
HermCMatrixT<REAL> &HermCMatrixT<REAL>::transpose() {
  return conj();
}

template<typename REAL>
const HermCMatrixT<REAL> &HermCMatrixT<REAL>::adjoint() const {
  return *this;
}

template<typename REAL>
HermCMatrixT<REAL> &HermCMatrixT<REAL>::Identity() {
  for (int i = 0; i < dim; i++)
    for (int j = i; j < dim; j++)
      (*this)(COMPLEX(i == j), i, j);
  return *this;
}

template<>
HermCMatrixT<double> &HermCMatrixT<double>::Invert() {
  std::complex<double> A[dim * dim];
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      A[i * dim + j] = (*this)(i, j);

  int I[dim];
  int info;
  int lwork = 2 * dim;
  std::complex<double> work[2 * lwork];
  char uplo = 'L';
  zsytrf_(&uplo, &dim, A, &dim, I, work, &lwork, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }
  zsytri_(&uplo, &dim, A, &dim, I, work, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }

  for (int i = 0; i < dim; ++i)
    for (int j = i; j < dim; ++j)
      (*this)(A[i * dim + j], i, j);

  return *this;
}

template<typename REAL>
Saver &HermCMatrixT<REAL>::save(Saver &saver) const {
  saver << Dim();
  for (int i = 0; i < size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &HermCMatrixT<REAL>::load(Loader &loader) {
  int N;
  loader >> N;
  resize(N);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template class HermCMatrixT<double>;

#ifdef BUILD_IA

template class HermCMatrixT<IAInterval>;

HermCMatrix mid(const IAHermCMatrix &IA_A) {
  HermCMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A(mid(IA_A(i, j)), i, j);
  return A;
}

#endif // BUILD_IA