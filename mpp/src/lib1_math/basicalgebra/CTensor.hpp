#ifndef MLUQ_CTENSOR_H
#define MLUQ_CTENSOR_H

#include "CMatrix.hpp"
#include "RTensor.hpp"

template<typename REAL = double>
class CTensorT {
  using COMPLEX = COMPLEX_TYPE<REAL>;
protected:
  std::vector<COMPLEX> z{};
  int NFirstDimension;
  int NSecondDimension;
  int NThirdDimension;
public:
  CTensorT() : z(0, COMPLEX{}), NFirstDimension(0), NSecondDimension(0), NThirdDimension(0) {}

  explicit constexpr CTensorT(int dim) :
      z(dim * dim * dim, COMPLEX{}), NFirstDimension(dim), NSecondDimension(dim),
      NThirdDimension(dim) {}

  explicit constexpr CTensorT(const COMPLEX &b, int dim) :
      z(dim * dim * dim, b), NFirstDimension(dim), NSecondDimension(dim), NThirdDimension(dim) {}

  constexpr CTensorT(int FirstDimension, int SecondDimension, int ThirdDimension) :
      z(FirstDimension * SecondDimension * ThirdDimension, COMPLEX{}),
      NFirstDimension(FirstDimension), NSecondDimension(SecondDimension),
      NThirdDimension(ThirdDimension) {}

  constexpr CTensorT(const COMPLEX &b, int FirstDimension, int SecondDimension,
                     int ThirdDimension) :
      z(FirstDimension * SecondDimension * ThirdDimension, b), NFirstDimension(FirstDimension),
      NSecondDimension(SecondDimension), NThirdDimension(ThirdDimension) {}

  template<typename COMPLEX1>
  CTensorT(const std::initializer_list<std::initializer_list<std::initializer_list<COMPLEX1>>> &v) {
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

  template<typename COMPLEX1>
  CTensorT(const std::vector<std::vector<std::vector<COMPLEX1>>> &v) {
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
  explicit CTensorT(const CTensorT<REAL1> &A) :
      z(A.begin(), A.end()), NFirstDimension(A.FirstDimension()),
      NSecondDimension(A.SecondDimension()), NThirdDimension(A.ThirdDimension()) {}

  template<typename REAL1>
  explicit CTensorT(const RTensorT<REAL1> &A) :
      z(A.begin(), A.end()), NFirstDimension(A.FirstDimension()),
      NSecondDimension(A.SecondDimension()), NThirdDimension(A.ThirdDimension()) {}

  template<typename REAL1>
  explicit CTensorT(const CTensorT<REAL1> &A, const CTensorT<REAL1> &B) :
      z(A.begin(), A.end()), NFirstDimension(A.FirstDimension()),
      NSecondDimension(A.SecondDimension()), NThirdDimension(A.ThirdDimension()) {
    if (NFirstDimension != B.FirstDimension() || NSecondDimension != B.SecondDimension()
        || NThirdDimension != B.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] *= B.Data()[i];
  }

  template<typename REAL1>
  explicit CTensorT(const RTensorT<REAL1> &A, const CTensorT<REAL1> &B) :
      z(B.begin(), B.end()), NFirstDimension(B.FirstDimension()),
      NSecondDimension(B.SecondDimension()), NThirdDimension(B.ThirdDimension()) {
    if (NFirstDimension != A.FirstDimension() || NSecondDimension != A.SecondDimension()
        || NThirdDimension != A.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] *= A.Data()[i];
  }

  template<typename REAL1>
  explicit CTensorT(const CTensorT<REAL1> &A, const RTensorT<REAL1> &B) :
      z(A.begin(), A.end()), NFirstDimension(A.FirstDimension()),
      NSecondDimension(A.SecondDimension()), NThirdDimension(A.ThirdDimension()) {
    if (NFirstDimension != B.FirstDimension() || NSecondDimension != B.SecondDimension()
        || NThirdDimension != B.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] *= B.Data()[i];
  }

  template<typename REAL1>
  CTensorT &operator=(const CTensorT<REAL1> &A) {
    NFirstDimension = A.FirstDimension();
    NSecondDimension = A.SecondDimension();
    NThirdDimension = A.ThirdDimension();

    z = std::vector<COMPLEX>(A.begin(), A.end());
    return *this;
  }

  CTensorT &operator=(const COMPLEX &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] = b;
    return *this;
  }

  template<typename REAL1>
  CTensorT &operator+=(const CTensorT<REAL1> &A) {
    if (NFirstDimension != A.FirstDimension() || NSecondDimension != A.SecondDimension()
        || NThirdDimension != A.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] += A.Data()[i];
    return *this;
  }

  template<typename REAL1>
  CTensorT &operator-=(const CTensorT<REAL1> &A) {
    if (NFirstDimension != A.FirstDimension() || NSecondDimension != A.SecondDimension()
        || NThirdDimension != A.ThirdDimension())
      THROW("Size does not fit!")
    for (int i = 0; i < z.size(); ++i)
      z[i] -= A.Data()[i];
    return *this;
  }

  CTensorT &operator*=(const COMPLEX &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] *= b;
    return *this;
  }

  CTensorT &operator/=(const COMPLEX &b) {
    for (int i = 0; i < z.size(); ++i)
      z[i] /= b;
    return *this;
  }

  // Todo: Next to operators are only a draft

  auto operator[](int i) { return z.begin() + i * NFirstDimension; } //????

  auto operator[](int i) const { return z.begin() + i * NFirstDimension; } //????

  COMPLEX &operator()(int FirstComponent, int SecondComponent, int ThirdComponent) {
    return z[FirstComponent + SecondComponent * NFirstDimension
             + ThirdComponent * NFirstDimension * NSecondDimension];
  }

  // Todo: Arguments need to be refactored. Is it useful to have these functions?

  CMatrixT<REAL> FirstComponent(int i) const {
    CMatrixT<REAL> c(NSecondDimension, NThirdDimension);
    for (int j = 0; j < NSecondDimension; ++j) {
      for (int k = 0; k < NThirdDimension; ++k) {
        c(j, k) = z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension];
      }
    }
    return c;
  }

  CMatrixT<REAL> SecondComponent(int j) const {
    CMatrixT<REAL> c(NFirstDimension, NThirdDimension);
    for (int i = 0; i < NFirstDimension; ++i) {
      for (int k = 0; k < NThirdDimension; ++k) {
        c(i, k) = z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension];
      }
    }
    return c;
  }

  CMatrixT<REAL> ThirdComponent(int k) const {
    CMatrixT<REAL> c(NFirstDimension, NSecondDimension);
    for (int i = 0; i < NFirstDimension; ++i) {
      for (int j = 0; j < NSecondDimension; ++j) {
        c(i, j) = z[i + j * NFirstDimension + k * NFirstDimension * NSecondDimension];
      }
    }
    return c;
  }

  template<typename REAL1>
  void InsertMatrixThirdDimension(const CMatrixT<REAL1> &v, int ThirdComponent, int startFirst = 0,
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

  const std::vector<COMPLEX> &Data() const { return z; }

  std::vector<COMPLEX> &Data() { return z; }

  RTensorT<REAL> real() const;

  RTensorT<REAL> imag() const;

  COMPLEX Mean() const;

  REAL Variance() const;
};

template<typename REAL>
inline CTensorT<REAL> operator-(const CTensorT<REAL> &a) {
  CTensorT<REAL> c(a);
  return c *= -1;
}

template<typename REAL>
inline CTensorT<REAL> operator+(const CTensorT<REAL> &a, const CTensorT<REAL> &b) {
  CTensorT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CTensorT<REAL> operator+(const CTensorT<REAL> &a, const RTensorT<REAL> &b) {
  CTensorT<REAL> c(b);
  return c += a;
}

template<typename REAL>
inline CTensorT<REAL> operator+(const RTensorT<REAL> &a, const CTensorT<REAL> &b) {
  CTensorT<REAL> c(a);
  return c += b;
}

template<typename REAL>
inline CTensorT<REAL> operator-(const CTensorT<REAL> &a, const CTensorT<REAL> &b) {
  CTensorT<REAL> c(a);
  return c -= b;
}

template<typename REAL>
inline CTensorT<REAL> operator-(const CTensorT<REAL> &a, const RTensorT<REAL> &b) {
  CTensorT<REAL> c(a);
  CTensorT<REAL> d(b);
  return c -= d;
}

template<typename REAL>
inline CTensorT<REAL> operator-(const RTensorT<REAL> &a, const CTensorT<REAL> &b) {
  CTensorT<REAL> c(a);
  return c -= b;
}

template<typename REAL> // elementwise multiplication
inline CTensorT<REAL> operator*(const CTensorT<REAL> &a, const CTensorT<REAL> &b) {
  return CTensorT<REAL>(a, b);
}

template<typename REAL> // elementwise multiplication
inline CTensorT<REAL> operator*(const CTensorT<REAL> &a, const RTensorT<REAL> &b) {
  return CTensorT<REAL>(a, b);
}

template<typename REAL> // elementwise multiplication
inline CTensorT<REAL> operator*(const RTensorT<REAL> &a, const CTensorT<REAL> &b) {
  return CTensorT<REAL>(a, b);
}

template<typename REAL, typename COMPLEX1>
inline CTensorT<REAL> operator*(const CTensorT<REAL> &a, const COMPLEX1 &b) {
  CTensorT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename COMPLEX1>
inline CTensorT<REAL> operator*(const COMPLEX1 &b, const CTensorT<REAL> &a) {
  CTensorT<REAL> c(a);
  return c *= b;
}

template<typename REAL, typename COMPLEX1>
inline CTensorT<REAL> operator/(const CTensorT<REAL> &a, const COMPLEX1 &b) {
  CTensorT<REAL> c(a);
  return c /= b;
}

template<typename REAL>
bool operator==(const CTensorT<REAL> &A, const CTensorT<REAL> &B) {
  if (A.FirstDimension() != B.FirstDimension() || A.SecondDimension() != B.SecondDimension()
      || A.ThirdDimension() != B.ThirdDimension())
    return false;
  return (A.real() == B.real()) && (A.imag() == B.imag());
}

template<typename REAL>
bool operator==(const CTensorT<REAL> &A, const RTensorT<REAL> &B) {
  return A == CTensorT<REAL>(B);
}

template<typename REAL>
bool operator==(const RTensorT<REAL> &A, const CTensorT<REAL> &B) {
  return CTensorT<REAL>(A) == B;
}

using CTensor = CTensorT<double>;

#endif // MLUQ_CTENSOR_HPP
