#include "CMatrix.hpp"
#include "CFunctions.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

#include <cmath>

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
// LU factorization
extern "C" void zgetrf_(int *M, int *N, void *A, int *LDA, int *IPIV, int *INFO);
// Invert matrix using LU factorization
extern "C" void zgetri_(int *N, void *A, int *LDA, int *IPIV, void *WORK, int *LWORK, int *INFO);

//==================================================================================================

template<typename REAL>
CVectorT<REAL> CMatrixT<REAL>::diag() {
  CVectorT<REAL> d(std::min(Nrows, Ncols));
  for (int i = 0; i < std::min(Nrows, Ncols); ++i)
    d[i] = z[i * Ncols + i];
  return d;
}

template<typename REAL>
void CMatrixT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
CMatrixT<REAL> &CMatrixT<REAL>::conj() {
  for (int i = 0; i < z.size(); ++i)
    z[i] = baConj(z[i]);
  return *this;
}

template<typename REAL>
RMatrixT<REAL> CMatrixT<REAL>::real() const {
  RMatrixT<REAL> re(Nrows, Ncols);
  for (int i = 0; i < z.size(); ++i)
    re.Data()[i] = baReal(z[i]);
  return re;
}

template<typename REAL>
RMatrixT<REAL> CMatrixT<REAL>::imag() const {
  RMatrixT<REAL> im(Nrows, Ncols);
  for (int i = 0; i < z.size(); ++i)
    im.Data()[i] = baImag(z[i]);
  return im;
}

template<typename REAL>
CMatrixT<REAL> &CMatrixT<REAL>::transpose() {
  CMatrixT<REAL> b(*this);
  resize(b.Ncols, b.Nrows);
  for (int i = 0; i < Nrows; ++i)
    for (int j = 0; j < Ncols; ++j)
      z[i * Ncols + j] = b[j][i];
  return *this;
}

template<typename REAL>
COMPLEX_TYPE<REAL> CMatrixT<REAL>::Mean() const {
  COMPLEX mean{};
  for (auto const &value_vec : z)
    mean += value_vec;
  return mean / z.size();
}

template<typename REAL>
REAL CMatrixT<REAL>::Variance() const {
  REAL var_real = real().Variance();
  REAL var_imag = imag().Variance();
  return var_real + var_imag;
}

template<typename REAL>
CMatrixT<REAL> &CMatrixT<REAL>::adjoint() {
  transpose();
  conj();
  return *this;
}

template<typename REAL>
CMatrixT<REAL> &CMatrixT<REAL>::Identity() {
  for (int i = 0; i < Nrows; i++)
    for (int j = 0; j < Ncols; j++)
      z[i * Ncols + j] = COMPLEX(i == j);
  return *this;
}

template<>
CMatrixT<double> &CMatrixT<double>::Invert() {
  if (Nrows != Ncols) THROW("Nonquadratic matix in invert")

  int I[Nrows];
  int info;
  int lwork = 2 * Nrows;
  std::complex<double> work[2 * lwork];
  zgetrf_(&Nrows, &Nrows, z.data(), &Nrows, I, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }
  zgetri_(&Nrows, z.data(), &Nrows, I, work, &lwork, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }

  return *this;
}

template<typename REAL>
Saver &CMatrixT<REAL>::save(Saver &saver) const {
  saver << Nrows << Ncols;
  for (int i = 0; i < z.size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &CMatrixT<REAL>::load(Loader &loader) {
  int r, c;
  loader >> r >> c;
  resize(r, c);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template<typename REAL>
void CMatrixT<REAL>::resize(int rows, int cols) {
  z = std::vector<COMPLEX>(rows * cols, COMPLEX{});
  Nrows = rows;
  Ncols = cols;
}

template class CMatrixT<double>;

#ifdef BUILD_IA

template class CMatrixT<IAInterval>;

CMatrix mid(const IACMatrix &IA_A) {
  CMatrix A(IA_A.rows(), IA_A.cols());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A[i][j] = mid(IA_A[i][j]);
  return A;
}

#endif // BUILD_IA
