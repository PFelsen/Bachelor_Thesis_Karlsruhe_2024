#ifndef MLUQ_RTENSOR_H
#define MLUQ_RTENSOR_H

#include "RMatrix.hpp"

template<typename REAL = double>
class RTensorT {
  std::vector<REAL> z{};
  int NFirstDimension;
  int NSecondDimension;
  int NThirdDimension;
public:
  RTensorT() : z(0, REAL{}), NFirstDimension(0), NSecondDimension(0), NThirdDimension(0) {}

  explicit constexpr RTensorT(int dim) :
      z(dim * dim * dim, REAL{}), NFirstDimension(dim), NSecondDimension(dim),
      NThirdDimension(dim) {}

  explicit constexpr RTensorT(const REAL &b, int dim) :
      z(dim * dim * dim, b), NFirstDimension(dim), NSecondDimension(dim), NThirdDimension(dim) {}

  constexpr RTensorT(int FirstDimension, int SecondDimension, int ThirdDimension) :
      z(FirstDimension * SecondDimension * ThirdDimension, REAL{}), NFirstDimension(FirstDimension),
      NSecondDimension(SecondDimension), NThirdDimension(ThirdDimension) {}

  constexpr RTensorT(const REAL &b, int FirstDimension, int SecondDimension, int ThirdDimension) :
      z(FirstDimension * SecondDimension * ThirdDimension, b), NFirstDimension(FirstDimension),
      NSecondDimension(SecondDimension), NThirdDimension(ThirdDimension) {}

  template<typename REAL1>
  RTensorT(const std::initializer_list<std::initializer_list<std::initializer_list<REAL1>>> &v) {
    NFirstDimension = (int)(v.begin())->begin()->size();
    NSecondDimension = (int)(v.begin())->size();
    NThirdDimension = (int)v.size();
    z.resize(NFirstDimension * NSecondDimension * NThirdDimension);
    for (int i = 0; i < NFirstDimension; i++) {
      for (int j = 0; j < NSecondDimension; j++) {
        for (int k = 0; k < NThirdDimension; k++) {
          z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension] =
              (((v.begin() + i)->begin() + j)->begin())[k];
        }
      }
    }
  }

  template<typename REAL1>
  RTensorT(const std::vector<std::vector<std::vector<REAL1>>> &v) {
    NFirstDimension = (int)(v.begin())->begin()->size();
    NSecondDimension = (int)(v.begin())->size();
    NThirdDimension = (int)v.size();
    z.resize(NFirstDimension * NSecondDimension * NThirdDimension);
    for (int i = 0; i < NFirstDimension; i++) {
      for (int j = 0; j < NSecondDimension; j++) {
        for (int k = 0; k < NThirdDimension; k++) {
          z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension] =
              (((v.begin() + i)->begin() + j)->begin())[k];
        }
      }
    }
  }

  template<typename REAL1>
  explicit RTensorT(const RTensorT<REAL1> &A) :
      z(A.begin(), A.end()), NFirstDimension(A.FirstDimension()),
      NSecondDimension(A.SecondDimension()), NThirdDimension(A.ThirdDimension()) {}

  template<typename REAL1>
  explicit RTensorT(const RTensorT<REAL1> &A, const RTensorT<REAL1> &B) :
      z(A.begin(), A.end()), NFirstDimension(A.FirstDimension()),
      NSecondDimension(A.SecondDimension()), NThirdDimension(A.ThirdDimension()) {
    if (NFirstDimension != B.FirstDimension() || NSecondDimension != B.SecondDimension()
        || NThirdDimension != B.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] *= B.Data()[i];
  }

  template<typename REAL1>
  RTensorT &operator=(const RTensorT<REAL1> &A) {
    NFirstDimension = A.FirstDimension();
    NSecondDimension = A.SecondDimension();
    NThirdDimension = A.ThirdDimension();

    z = std::vector<REAL>(A.begin(), A.end());
    return *this;
  }

  RTensorT &operator=(const REAL &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = b;
    return *this;
  }

  template<typename REAL1>
  RTensorT &operator+=(const RTensorT<REAL1> &A) {
    if (NFirstDimension != A.FirstDimension() || NSecondDimension != A.SecondDimension()
        || NThirdDimension != A.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A.Data()[i];
    return *this;
  }

  template<typename REAL1>
  RTensorT &operator-=(const RTensorT<REAL1> &A) {
    if (NFirstDimension != A.FirstDimension() || NSecondDimension != A.SecondDimension()
        || NThirdDimension != A.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A.Data()[i];
    return *this;
  }

  RTensorT &operator*=(const REAL &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= b;
    return *this;
  }

  RTensorT &operator/=(const REAL &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= b;
    return *this;
  }

  // Todo: Next to operators are only a draft

  auto operator[](int i) { return z.begin() + i * NFirstDimension; } //????

  auto operator[](int i) const { return z.begin() + i * NFirstDimension; } //????

  REAL &operator()(int FirstComponent, int SecondComponent, int ThirdComponent) {
    return z[FirstComponent + SecondComponent * NFirstDimension
             + ThirdComponent * NFirstDimension * NSecondDimension];
  }

  auto &operator()(int FirstComponent, int SecondComponent, int ThirdComponent) const {
    return z[FirstComponent + SecondComponent * NFirstDimension
             + ThirdComponent * NFirstDimension * NSecondDimension];
  }

  // Todo: Arguments need to be refactored. Is it useful to have these functions?

  RMatrixT<REAL> FirstComponent(int i) const {
    RMatrixT<REAL> c(NSecondDimension, NThirdDimension);
    for (int j = 0; j < NSecondDimension; ++j) {
      for (int k = 0; k < NThirdDimension; ++k) {
        c(j, k) = z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension];
      }
    }
    return c;
  }

  RMatrixT<REAL> SecondComponent(int j) const {
    RMatrixT<REAL> c(NFirstDimension, NThirdDimension);
    for (int i = 0; i < NFirstDimension; ++i) {
      for (int k = 0; k < NThirdDimension; ++k) {
        c(i, k) = z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension];
      }
    }
    return c;
  }

  RMatrixT<REAL> ThirdComponent(int k) const {
    RMatrixT<REAL> c(NFirstDimension, NSecondDimension);
    for (int i = 0; i < NFirstDimension; ++i) {
      for (int j = 0; j < NSecondDimension; ++j) {
        c(i, j) = z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension];
      }
    }
    return c;
  }

  template<typename REAL1>
  void InsertMatrixThirdDimension(const RMatrixT<REAL1> &v, int ThirdComponent, int startFirst = 0,
                                  int startSecond = 0) {
    if (ThirdComponent >= NThirdDimension || startFirst + v.cols() > NFirstDimension
        || startSecond + v.rows() > NSecondDimension) {
      THROW("Cannot insert Matrix: size of v too large")
    }
    for (int i = 0; i < v.cols(); ++i) {
      for (int j = 0; j < v.rows(); ++j) {
        z[i + startFirst + (j + startSecond) * NFirstDimension
          + ThirdComponent * NFirstDimension * NSecondDimension] = v(i, j);
      }
    }
  }

  int size() const { return z.size(); }

  int FirstDimension() const { return NFirstDimension; }

  int SecondDimension() const { return NSecondDimension; }

  int ThirdDimension() const { return NThirdDimension; }

  auto begin() const { return z.begin(); }

  auto begin() { return z.begin(); }

  auto end() const { return z.end(); }

  auto end() { return z.end(); }

  void resize(int FirstComponents, int SecondComponents, int ThirdComponents);

  void resize(int dim) { resize(dim, dim, dim); }

  const std::vector<REAL> &Data() const { return z; }

  std::vector<REAL> &Data() { return z; }

  REAL Mean() const;

  REAL Variance() const;
};

template<typename REAL>
inline RTensorT<REAL> operator-(const RTensorT<REAL> &a) {
  RTensorT<REAL> c(a);
  return c *= -1;
}

template<typename REAL>
inline RTensorT<REAL> operator+(const RTensorT<REAL> &a, const RTensorT<REAL> &b) {
  RTensorT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline RTensorT<REAL> operator-(const RTensorT<REAL> &a, const RTensorT<REAL> &b) {
  RTensorT<REAL> c(a);
  return c -= b;
}

template<typename REAL> // elementwise multiplication
inline RTensorT<REAL> operator*(const RTensorT<REAL> &a, const RTensorT<REAL> &b) {
  return RTensorT<REAL>(a, b);
}

template<typename REAL, typename REAL1>
inline RTensorT<REAL> operator*(const RTensorT<REAL> &a, const REAL1 &b) {
  RTensorT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename REAL1>
inline RTensorT<REAL> operator*(const REAL1 &b, const RTensorT<REAL> &a) {
  RTensorT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename REAL1>
inline RTensorT<REAL> operator/(const RTensorT<REAL> &a, const REAL1 &b) {
  RTensorT<REAL> c(a);
  return c /= b;
}

template<typename REAL>
bool operator==(const RTensorT<REAL> &A, const RTensorT<REAL> &B) {
  if (A.FirstDimension() != B.FirstDimension() || A.SecondDimension() != B.SecondDimension()
      || A.ThirdDimension() != B.ThirdDimension())
    return false;
  for (int i = 0; i < A.Data().size(); ++i)
    if (!mpp_ba::isNear(A.Data()[i], B.Data()[i])) return false;
  return true;
}

using RTensor = RTensorT<double>;

#endif // MLUQ_RTENSOR_H
