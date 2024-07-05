#include "AntisymRMatrix.hpp"
#include "Parallel.hpp"
#include "RMatrix.hpp"
#include "SaveLoad.hpp"

#include <cmath>

template<typename REAL>
AntisymRMatrixT<REAL>::AntisymRMatrixT(int dim) : dim(dim), z((dim * (dim - 1)) / 2, REAL{}) {}

template<typename REAL>
REAL AntisymRMatrixT<REAL>::operator()(int i, int j) {
  if (i == j) return REAL{};
  if (i > j) return z[(i * (i - 1)) / 2 + j];
  else return -z[(j * (j - 1)) / 2 + i];
}

template<typename REAL>
REAL AntisymRMatrixT<REAL>::operator()(int i, int j) const {
  if (i == j) return REAL{};
  if (i > j) return z[(i * (i - 1)) / 2 + j];
  else return -z[(j * (j - 1)) / 2 + i];
}

template<typename REAL>
void AntisymRMatrixT<REAL>::resize(int dim) {
  z.resize((dim * (dim - 1)) / 2);
  this->dim = dim;
}

template<typename REAL>
AntisymRMatrixT<REAL> &AntisymRMatrixT<REAL>::transpose() {
  for (int i = 0; i < z.size(); ++i)
    z[i] *= -1.0;
  return *this;
}

template<>
AntisymRMatrixT<double> &AntisymRMatrixT<double>::Invert() {
  if (dim % 2 == 1)
    Exit("Antisymmetric matrix not invertible for dim= " + std::to_string(dim))

        double a[dim * dim];
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      a[i * dim + j] = (*this)(i, j);

  int I[dim];
  int info;
  int lwork = 2 * dim;
  double work[2 * lwork];
  dgetrf_(&dim, &dim, a, &dim, I, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }
  dgetri_(&dim, a, &dim, I, work, &lwork, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }

  for (int i = 1; i < dim; ++i)
    for (int j = 0; j < i; ++j)
      (*this)(a[i * dim + j], i, j);

  return *this;
}

template<typename REAL>
void AntisymRMatrixT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
Saver &AntisymRMatrixT<REAL>::save(Saver &saver) const {
  saver << Dim();
  for (int i = 0; i < z.size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &AntisymRMatrixT<REAL>::load(Loader &loader) {
  int N;
  loader >> N;
  resize(N);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template class AntisymRMatrixT<double>;

#ifdef BUILD_IA

template class AntisymRMatrixT<IAInterval>;

AntisymRMatrix mid(const IAAntisymRMatrix &IA_A) {
  AntisymRMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A(mid(IA_A(i, j)), i, j);
  return A;
}

RMatrix sup(const IAAntisymRMatrix &IA_A) {
  RMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A[i][j] = sup(IA_A(i, j));
  return A;
}

RMatrix inf(const IAAntisymRMatrix &IA_A) {
  RMatrix A(IA_A.Dim());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A[i][j] = inf(IA_A(i, j));
  return A;
}

#endif // BUILD_IA