#ifndef QUADRATURESYM_HPP
#define QUADRATURESYM_HPP

#include "QuadratureT.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QemptyT : public QuadratureT<T, sDim, tDim> {
public:
  QemptyT() : QuadratureT<T, sDim, tDim>("Qempty") {}
};

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QpointT : public QuadratureT<T, sDim, tDim> {
public:
  QpointT();
};

typedef QpointT<> Qpoint;

#ifdef BUILD_IA

typedef QpointT<IAInterval, SpaceDimension, TimeDimension> IAQpoint;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QintT : public QuadratureT<T, sDim, tDim> {
  void initDataPair(int cnt, T x, T w) {
    this->z[cnt] = PointT<T, sDim, tDim>(x);
    this->w[cnt] = w;
  }

  void setName(int exactUpTo);
public:
  QintT(int exactUpTo);
};

typedef QintT<> Qint;

#ifdef BUILD_IA

typedef QintT<IAInterval, SpaceDimension, TimeDimension> IAQint;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QtriSymT : public QuadratureT<T, sDim, tDim> {
  void initDataPair(int cnt, T x, T y, T w) {
    this->z[cnt] = PointT<T, sDim, tDim>(x, y);
    this->w[cnt] = w;
  }

  void setName(int exactUpTo);
public:
  QtriSymT(int exactUpTo);
};

typedef QtriSymT<> QtriSym;

#ifdef BUILD_IA

typedef QtriSymT<IAInterval, SpaceDimension, TimeDimension> IAQtriSym;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QquadT : public QTensorT<T, sDim, tDim> {
private:
  std::vector<QuadratureT<T, sDim, tDim>> createQuads(int &exactUpTo, int &exactUpToY) {
    if (exactUpToY < 0) exactUpToY = exactUpTo;
    return {QintT<T, sDim, tDim>(exactUpTo), QintT<T, sDim, tDim>(exactUpToY)};
  }
public:
  QquadT(int exactUpTo, int exactUpToY);
};

typedef QquadT<> Qquad;

#ifdef BUILD_IA

typedef QquadT<IAInterval, SpaceDimension, TimeDimension> IAQquad;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QtetSymT : public QuadratureT<T, sDim, tDim> {
  void initDataPair(int cnt, T x, T y, T z, T w) {
    this->z[cnt] = PointT<T, sDim, tDim>(x, y, z);
    this->w[cnt] = w;
  }

  void setName(int exactUpTo);
public:
  QtetSymT(int exactUpTo);
};

typedef QtetSymT<> QtetSym;

#ifdef BUILD_IA

typedef QtetSymT<IAInterval, SpaceDimension, TimeDimension> IAQtetSym;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QpriSymT : public QTensorT<T, sDim, tDim> {
private:
  std::vector<QuadratureT<T, sDim, tDim>> createQuads(int &exactUpTo, int &exactUpToY) {
    if (exactUpToY < 0) exactUpToY = exactUpTo;
    return {QtriSymT<T, sDim, tDim>(exactUpTo), QintT<T, sDim, tDim>(exactUpToY)};
  }
public:
  QpriSymT(int exactUpToTriangle, int exactUpToZ);
};

typedef QpriSymT<> QpriSym;

#ifdef BUILD_IA

typedef QpriSymT<IAInterval, SpaceDimension, TimeDimension> IAQpriSym;

#endif

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QhexT : public QTensorT<T, sDim, tDim> {
private:
  std::vector<QuadratureT<T, sDim, tDim>> createQuads(int &exactUpTo, int &exactUpToY,
                                                      int &exactUpToZ) {
    if (exactUpToY < 0 || exactUpToZ < 0) exactUpToY = exactUpToZ = exactUpTo;
    return {QintT<T, sDim, tDim>(exactUpTo), QintT<T, sDim, tDim>(exactUpToY),
            QintT<T, sDim, tDim>(exactUpToZ)};
  }
public:
  QhexT(int exactUpTo, int exactUpToY, int exactUpToZ);
};

typedef QhexT<double, SpaceDimension, TimeDimension> Qhex;

#ifdef BUILD_IA

typedef QhexT<IAInterval, SpaceDimension, TimeDimension> IAQhex;

#endif

#endif // QUADRATURESYM_HPP
