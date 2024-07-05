#include "RMatrix.hpp"
#include "GlobalDefinitions.hpp"
#include "Lapdef.hpp"
#include "Parallel.hpp"
#include "SaveLoad.hpp"

#include <cmath>

//==================================================================================================
// BLAS/LAPACK routines
//==================================================================================================
// LU factorization
extern "C" void dgetrf_(int *M, int *N, void *A, int *LDA, int *IPIV, int *INFO);
// Invert matrix using LU factorization
extern "C" void dgetri_(int *N, void *A, int *LDA, int *IPIV, void *WORK, int *LWORK, int *INFO);
//==================================================================================================
extern "C" void dgesv_(int *N, int *Nrhs, void *A, int *LDA, int *IPIV, void *B, int *LDB,
                       int *INFO, int ch1, int ch2);

template<typename REAL>
REAL RMatrixT<REAL>::normSqr() const {
  REAL sum{};
  for (int i = 0; i < z.size(); ++i)
    sum += z[i] * z[i];
  return sum;
}

template<typename REAL>
REAL RMatrixT<REAL>::norm() const {
  return sqrt(normSqr());
}

template<typename REAL>
RVectorT<REAL> RMatrixT<REAL>::diag() {
  RVectorT<REAL> d(std::min(Nrows, Ncols));
  for (int i = 0; i < std::min(Nrows, Ncols); ++i)
    d[i] = z[i * Ncols + i];
  return d;
}

template<typename REAL>
void RMatrixT<REAL>::resize(int rows, int cols) {
  z = std::vector<REAL>(rows * cols, REAL{});
  Nrows = rows;
  Ncols = cols;
}

template<typename REAL>
void RMatrixT<REAL>::Accumulate(int commSplit) {
  PPM->SumOnCommSplit(z, commSplit);
}

template<typename REAL>
RMatrixT<REAL> &RMatrixT<REAL>::transpose() {
  RMatrixT<REAL> b(*this);
  this->resize(b.Ncols, b.Nrows);
  for (int i = 0; i < Ncols; ++i)
    for (int j = 0; j < Nrows; ++j)
      (*this)[j][i] = b[i][j];
  return *this;
}

template<typename REAL>
RMatrixT<REAL> &RMatrixT<REAL>::Identity() {
  for (int i = 0; i < Nrows; i++)
    for (int j = 0; j < Ncols; j++)
      z[i * Ncols + j] = REAL(i == j);
  return *this;
}

template<>
RMatrixT<double> &RMatrixT<double>::Invert() {
  if (Nrows != Ncols) THROW("Nonquadratic matix in invert")

  int I[Nrows];
  int info;
  int lwork = 2 * Nrows;
  double work[2 * lwork];
  dgetrf_(&Nrows, &Nrows, z.data(), &Nrows, I, &info);

  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }
  dgetri_(&Nrows, z.data(), &Nrows, I, work, &lwork, &info);
  if (info < 0) {
    THROW("Error in LAPACK: " + std::to_string(-info) + "-th argument had an illegal value")
  } else if (info > 0) {
    THROW("Error in LAPACK: A(" + std::to_string(-info) + "," + std::to_string(-info)
          + " is exactly zero.")
  }

  return *this;
}

void LUDecomp(int n, Scalar *Mat, int *ipv) {
  for (int i = 0; i < n; ++i) {
    ipv[i] = i;
  }
  for (int i = 0; i < n; ++i) {
    int k = i;
    double piv = abs(Mat[i * n + i]);
    for (int j = i + 1; j < n; ++j) {
      double sum = abs(Mat[j * n + i]);
      if (sum > piv) {
        k = j;
        piv = sum;
      }
    }
    if (k != i) {
      std::swap(ipv[i], ipv[k]);
      for (int j = 0; j < n; ++j) {
        std::swap(Mat[k * n + j], Mat[i * n + j]);
      }
    }
    Scalar dinv = Mat[i * n + i];
    if (abs(dinv) < 1e-12) THROW("Error in SmallMatrix Inversion")
    dinv = Mat[i * n + i] = 1.0 / dinv;
    for (int j = i + 1; j < n; ++j) {
      Scalar piv = (Mat[j * n + i] *= dinv);
      for (int k = i + 1; k < n; ++k)
        Mat[j * n + k] -= Mat[i * n + k] * piv;
    }
  }
}

void invertSmallMatrix(int n, Scalar *a) {
  int *ipv = new int[n];
  Scalar *Mat = new Scalar[n * n];
  Scalar *rhs = new Scalar[n];
  for (int i = 0; i < n * n; ++i) {
    Mat[i] = a[i];
  }
  LUDecomp(n, Mat, ipv);
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      rhs[i] = 0;
    }
    rhs[k] = 1.0;
    for (int i = 0; i < n; ++i) {
      Scalar sum = rhs[ipv[i]];
      for (int j = 0; j < i; ++j)
        sum -= Mat[i * n + j] * a[j * n + k];
      a[i * n + k] = sum; // Lii = 1
    }
    for (int i = n - 1; i >= 0; i--) {
      Scalar sum = a[i * n + k];
      for (int j = i + 1; j < n; ++j)
        sum -= Mat[i * n + j] * a[j * n + k];
      a[i * n + k] = sum * Mat[i * n + i]; // Uii = Inv(Mii)
    }
  }
  delete[] rhs;
  delete[] Mat;
  delete[] ipv;
}

void applySmallMatrix(int n, Scalar *u, const Scalar *a, const Scalar *b) {
  for (int i = 0; i < n; ++i) {
    Scalar s = 0;
    for (int k = 0; k < n; ++k)
      s += a[i * n + k] * b[k];
    u[i] = s;
  }
}

std::vector<double> getCoeff(int deg) {
  switch (deg) {
  case 3:
    return {120.0, 60.0, 12.0, 1.0};
  case 5:
    return {30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0};
  case 7:
    return {17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0};
  case 9:
    return {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0,
            2162160.0,     110880.0,     3960.0,       90.0,        1.0};
  default:
    THROW("Error" + std::to_string(deg));
  }
}

void GetPadeNomDenom13(RMatrix &Mat, RMatrix &A2) {
  RMatrix A4(A2 * A2);
  RMatrix A6(A2 * A4);
  std::vector<double> PadeCoeff13{64764752532480000.0,
                                  32382376266240000.0,
                                  7771770303897600.0,
                                  1187353796428800.0,
                                  129060195264000.0,
                                  10559470521600.0,
                                  670442572800.0,
                                  33522128640.0,
                                  1323241920.0,
                                  40840800.0,
                                  960960.0,
                                  16380.0,
                                  182.0,
                                  1.0};
  ;
  RMatrix U(PadeCoeff13[9] * A2);
  RMatrix V(PadeCoeff13[8] * A2);
  for (int i = 0; i < U.size(); ++i) {
    U.Data().data()[i] += PadeCoeff13[11] * A4.Data().data()[i];
    U.Data().data()[i] += PadeCoeff13[13] * A6.Data().data()[i];
    V.Data().data()[i] += PadeCoeff13[10] * A4.Data().data()[i];
    V.Data().data()[i] += PadeCoeff13[12] * A6.Data().data()[i];
  }
  U = A6 * U;
  V = A6 * V;
  for (int i = 0; i < U.size(); ++i) {
    U.Data().data()[i] += PadeCoeff13[7] * A6.Data().data()[i];
    U.Data().data()[i] += PadeCoeff13[5] * A4.Data().data()[i];
    U.Data().data()[i] += PadeCoeff13[3] * A2.Data().data()[i];
  }
  for (int col = 0; col < U.cols(); ++col) {
    U[col][col] += PadeCoeff13[1];
  }
  Mat = Mat * U;
  for (int i = 0; i < V.size(); ++i) {
    V.Data().data()[i] += PadeCoeff13[6] * A6.Data().data()[i];
    V.Data().data()[i] += PadeCoeff13[4] * A4.Data().data()[i];
    V.Data().data()[i] += PadeCoeff13[2] * A2.Data().data()[i];
  }
  for (int col = 0; col < V.cols(); ++col) {
    V[col][col] += PadeCoeff13[0];
  }
  A2 = V - Mat;
  Mat += V;
}

void GetPadeNomDenom(RMatrix &Mat, RMatrix &A2, int deg) {
  if (deg == 13) {
    GetPadeNomDenom13(Mat, A2);
    return;
  }
  auto Coeff = getCoeff(deg);
  RMatrix U(RVector(Coeff[1], Mat.rows()));
  RMatrix V(RVector(Coeff[0], Mat.rows()));

  for (int i = 0; i < U.size(); ++i) {
    U.Data().data()[i] += Coeff[3] * A2.Data().data()[i];
    V.Data().data()[i] += Coeff[2] * A2.Data().data()[i];
  }
  int m = (deg - 1) / 2;
  RMatrix MatTmp(A2);
  for (int i = 2; i <= m; ++i) {
    MatTmp = MatTmp * A2;
    for (int j = 0; j < U.size(); ++j) {
      V.Data().data()[j] += Coeff[2 * i] * MatTmp.Data().data()[j];
      U.Data().data()[j] += Coeff[2 * i + 1] * MatTmp.Data().data()[j];
    }
  }
  Mat = Mat * U;
  A2 = V - Mat;
  Mat += V;
}

void Pade(RMatrix &Mat, int deg) {
  RMatrix MatSquared(Mat * Mat);
  GetPadeNomDenom(Mat, MatSquared, deg);
  int n = Mat.rows();
  int *piv = new int[n];
  int info;
  dgetrf_(&n, &n, MatSquared.Data().data(), &n, piv, &info);
  const char *TRANS = "N";
  dgetrs_(TRANS, &n, &n, MatSquared.Data().data(), &n, piv, Mat.Data().data(), &n, &info);
}

template<>
RMatrixT<double> &RMatrixT<double>::Exp(bool deg13) {
  std::vector<double> ThetaMFull{
      1.495585217958292e-2, 2.539398330063230e-1, 9.504178996162932e-1,
      2.097847961257068e0,  5.371920351148152e0,
  };
  double mu = (*this).Trace() / Nrows;
  for (int i = 0; i < Nrows; ++i) {
    (*this)[i][i] -= mu;
  }
  //    double OneNorm = (*this).NormOne();
  double OneNorm = (*this).norm();
  std::vector<int> degs = {3, 5, 7, 9};
  if (deg13) degs.push_back(13);
  for (int i = 0; i < degs.size(); ++i) {
    if (OneNorm <= ThetaMFull[i]) {
      Pade(*this, degs[i]);
      (*this) *= exp(mu);
      return *this;
    }
  }
  int s = std::ceil(log2(OneNorm / ThetaMFull[degs.size() - 1]));
  (*this) *= (1.0 / pow(2.0, s));
  Pade(*this, degs[degs.size() - 1]);
  for (int i = 0; i < s; ++i)
    *this = (*this) * (*this);
  *this *= exp(mu);
  return *this;
}

void FwdBwdSubstitution(int n, double *Mat, int *piv, double *rhs) {
  for (int i = 0; i < n; ++i) {
    Scalar sum = rhs[piv[i]];
    for (int j = 0; j < i; ++j)
      sum -= Mat[i * n + j] * rhs[j];
    rhs[i] = sum; // Lii = 1
  }
  for (int i = n - 1; i >= 0; i--) {
    Scalar sum = rhs[i];
    for (int j = i + 1; j < n; ++j)
      sum -= Mat[i * n + j] * rhs[j];
    rhs[i] = sum * Mat[i * n + i]; // Uii = Inv(Mii)
  }
}

template<>
RMatrixT<double> &RMatrixT<double>::Exp(RVector &Evaluated, bool deg13) {
  std::vector<double> ThetaMFull{
      1.495585217958292e-2, 2.539398330063230e-1, 9.504178996162932e-1,
      2.097847961257068e0,  5.371920351148152e0,
  };
  double mu = (*this).Trace() / Nrows;
  for (int i = 0; i < Nrows; ++i) {
    (*this)[i][i] -= mu;
  }
  //    double OneNorm = (*this).NormOne();
  double OneNorm = (*this).norm();
  std::vector<int> degs = {3, 5, 7, 9};
  if (deg13) degs.push_back(13);
  for (int i = 0; i < degs.size(); ++i) {
    if (OneNorm <= ThetaMFull[i]) {
      RMatrix MatSquared((*this) * (*this));
      if (degs[i] == 13) {
        int n = (*this).rows();
        int *piv = new int[n];
        int info;
        GetPadeNomDenom13((*this), MatSquared);
        dgetrf_(&n, &n, MatSquared.Data().data(), &n, piv, &info);
        const char *TRANS = "T";
        int NRHS = 1;
        Evaluated = (*this) * Evaluated;
        dgetrs_(TRANS, &n, &NRHS, MatSquared.Data().data(), &n, piv, Evaluated.asVector().data(),
                &n, &info); // LAPACK Version
      } else {
        GetPadeNomDenom((*this), MatSquared, degs[i]);
        int n = (*this).rows();
        int *piv = new int[n];
        LUDecomp(n, MatSquared.Data().data(), piv);
        Evaluated = (*this) * Evaluated;
        FwdBwdSubstitution(n, MatSquared.Data().data(), piv, Evaluated.asVector().data());
      }
      Evaluated *= exp(mu);
      return *this;
    }
  }
  int s = std::ceil(log2(OneNorm / ThetaMFull[degs.size() - 1]));
  (*this) *= (1.0 / pow(2.0, s));
  RMatrix MatSquared((*this) * (*this));
  GetPadeNomDenom((*this), MatSquared, degs[degs.size() - 1]);
  int n = (*this).rows();
  int *piv2 = new int[n];
  int *piv = new int[n];

  int info;
  if (deg13) {
    dgetrf_(&n, &n, MatSquared.Data().data(), &n, piv, &info);
    int NRHS = 1;
    const char *TRANS = "T";
    for (int outer = 0; outer < pow(2, s); ++outer) {
      Evaluated = (*this) * Evaluated;
      dgetrs_(TRANS, &n, &NRHS, MatSquared.Data().data(), &n, piv, Evaluated.asVector().data(), &n,
              &info);
    }
  } else {
    LUDecomp(n, MatSquared.Data().data(), piv2);
    for (int outer = 0; outer < pow(2, s); ++outer) {
      Evaluated = (*this) * Evaluated;
      FwdBwdSubstitution(n, MatSquared.Data().data(), piv2, Evaluated.asVector().data());
    }
  }
  Evaluated *= exp(mu);
  return *this;
}

template<>
RMatrixT<double> &RMatrixT<double>::Phi1() {
  int K = (*this).rows();
  RMatrix XX(2 * K, 2 * K);
  for (int k = 0; k < K; ++k) {
    XX[k][K + k] = 1;
    for (int l = 0; l < K; ++l) {
      XX[K + k][K + l] = (*this)[k][l];
    }
  }
  XX.Exp();
  for (int k = 0; k < K; ++k)
    for (int l = 0; l < K; ++l)
      (*this)[k][l] = XX[k][K + l];
  return *this;
}

template<>
RMatrixT<double> &RMatrixT<double>::Phi1(RVector &Evaluated, bool deg13) {
  RMatrix H_k_k(0.0, (*this).cols() + 1);
  H_k_k.Insert((*this), 0, 0);
  for (int i = 0; i < (*this).cols(); ++i) {
    H_k_k[i][(*this).cols()] = Evaluated[i];
  }
  Evaluated = 0.0;
  Evaluated.push_back(1.0);
  H_k_k.Exp(Evaluated, deg13);
  Evaluated.removeLast();
  return *this;
}

template<>
RMatrixT<double> &RMatrixT<double>::Phi2() {
  int K = (*this).rows();
  RMatrix XX(2 * K, 2 * K);
  for (int k = 0; k < K; ++k) {
    XX[k][K + k] = 1;
    for (int l = 0; l < K; ++l) {
      XX[K + k][K + l] = (*this)[k][l];
    }
  }
  XX.Phi1();
  for (int k = 0; k < K; ++k)
    for (int l = 0; l < K; ++l)
      (*this)[k][l] = XX[k][K + l];
  return *this;
}

template<>
RMatrixT<double> &RMatrixT<double>::Phi3() {
  int K = (*this).rows();
  RMatrix XX(2 * K, 2 * K);
  for (int k = 0; k < K; ++k) {
    XX[k][K + k] = 1;
    for (int l = 0; l < K; ++l) {
      XX[K + k][K + l] = (*this)[k][l];
    }
  }
  XX.Phi2();
  for (int k = 0; k < K; ++k)
    for (int l = 0; l < K; ++l)
      (*this)[k][l] = XX[k][K + l];
  return *this;
}

template<typename REAL>
void RMatrixT<REAL>::SaddlePoint(const RMatrixT<REAL> &A, const RMatrixT<REAL> &B) {
  if (A.rows() == B.rows()) {
    resize(A.rows() + B.cols(), A.cols() + B.cols());
    for (int i = 0; i < A.rows(); ++i)
      for (int j = 0; j < A.cols(); ++j)
        (*this)[i][j] = A[i][j];
    for (int i = 0; i < B.cols(); ++i)
      for (int j = 0; j < B.rows(); ++j)
        (*this)[A.rows() + i][j] = (*this)[j][A.rows() + i] = B[j][i];
    for (int i = 0; i < B.cols(); ++i)
      for (int j = 0; j < B.cols(); ++j)
        (*this)[A.rows() + i][A.cols() + j] = 0;
  } else {
    resize(A.rows() + B.rows(), A.cols() + B.rows());
    for (int i = 0; i < A.rows(); ++i)
      for (int j = 0; j < A.cols(); ++j)
        (*this)[i][j] = A[i][j];
    for (int i = 0; i < B.rows(); ++i)
      for (int j = 0; j < B.cols(); ++j)
        (*this)[A.rows() + i][j] = (*this)[j][A.rows() + i] = B[i][j];
    for (int i = 0; i < B.rows(); ++i)
      for (int j = 0; j < B.rows(); ++j)
        (*this)[A.rows() + i][A.cols() + j] = 0;
  }
}

template<typename REAL>
Saver &RMatrixT<REAL>::save(Saver &saver) const {
  saver << Nrows << Ncols;
  for (int i = 0; i < z.size(); ++i)
    saver << z[i];
  return saver;
}

template<typename REAL>
Loader &RMatrixT<REAL>::load(Loader &loader) {
  int r, c;
  loader >> r >> c;
  resize(r, c);
  for (int i = 0; i < z.size(); ++i)
    loader >> z[i];
  return loader;
}

template<typename REAL>
REAL RMatrixT<REAL>::Mean() const {
  REAL mean{};
  for (auto const &value : z)
    mean += value;
  return mean / REAL(z.size());
}

template<typename REAL>
REAL RMatrixT<REAL>::Variance() const {
  REAL m2{};
  REAL diff{};
  REAL correction{};
  REAL mean = this->Mean();
  for (auto const &value : z) {
    diff = (value - mean);
    m2 += diff * diff;
    correction += diff;
  }
  return (m2 - correction * correction / z.size()) / (z.size() - 1);
}

template<typename REAL>
REAL RMatrixT<REAL>::Trace() const {
  REAL tr{};
  if (Nrows != Ncols) THROW("Trace only defined for square matrices!")
  for (int i = 0; i < Nrows; ++i) {
    tr += (*this)[i][i];
  }
  return tr;
}

template<typename REAL>
REAL RMatrixT<REAL>::NormOne() const {
  REAL norm{};
  REAL tmp{};
  for (int j = 0; j < Ncols; ++j) {
    for (int i = 0; i < Nrows; ++i) {
      tmp += abs((*this)[i][j]);
    }
    norm = max(norm, tmp);
  }
  return norm;
}

template<typename REAL>
REAL RMatrixT<REAL>::NormInfty() const {
  REAL norm{};
  REAL tmp{};
  for (int i = 0; i < Nrows; ++i) {
    for (int j = 0; j < Ncols; ++j) {
      tmp += abs((*this)[i][j]);
    }
    norm = max(norm, tmp);
  }
  return norm;
}

template class RMatrixT<double>;

#ifdef BUILD_IA

template class RMatrixT<IAInterval>;

RMatrix mid(const IARMatrix &IA_A) {
  RMatrix A(IA_A.rows(), IA_A.cols());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A[i][j] = mid(IA_A[i][j]);
  return A;
}

RMatrix sup(const IARMatrix &IA_A) {
  RMatrix A(IA_A.rows(), IA_A.cols());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A[i][j] = sup(IA_A[i][j]);
  return A;
}

RMatrix inf(const IARMatrix &IA_A) {
  RMatrix A(IA_A.rows(), IA_A.cols());
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.cols(); ++j)
      A[i][j] = inf(IA_A[i][j]);
  return A;
}

#endif // BUILD_IA