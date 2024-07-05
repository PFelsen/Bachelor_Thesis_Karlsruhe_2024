#include "QuadratureSym.hpp"
#include "QuadratureDataInterval.hpp"
#include "QuadratureDataTetrahedron.hpp"
#include "QuadratureDataTriangle.hpp"

template<>
QpointT<double, SpaceDimension, TimeDimension>::QpointT() : Quadrature("Qpoint", 1) {
  z[0] = Point(1.0, 0.0, 0.0);
  w[0] = 1.0;
}

#ifdef BUILD_IA

template<>
QpointT<IAInterval, SpaceDimension, TimeDimension>::QpointT() : IAQuadrature("IAQpoint", 1) {
  z[0] = IAPoint(IAInterval(1.0), IAInterval(0.0));
  w[0] = IAInterval(1.0);
}

#endif

template<>
void QintT<double, SpaceDimension, TimeDimension>::setName(int exactUpTo) {
  name = "Qint: exact up to degree " + std::to_string(exactUpTo);
}

template<typename T, int sDim, int tDim>
QintT<T, sDim, tDim>::QintT(int exactUpTo) {
  QintSymT_data<T> qData = GetQintSymT_data<T>(exactUpTo);
  this->resize(qData.Size());
  int cnt = 0;
  if (qData.d0.size() > 0) { initDataPair(cnt++, T(0.5), qData.d0[0].weight); }
  for (int i = 0; i < qData.d1.size(); ++i) {
    initDataPair(cnt++, qData.d1[i].A, qData.d1[i].weight);
    initDataPair(cnt++, qData.d1[i].B, qData.d1[i].weight);
  }
  setName(exactUpTo);
}

template class QintT<double, SpaceDimension, TimeDimension>;

#ifdef BUILD_IA

template<>
void QintT<IAInterval, SpaceDimension, TimeDimension>::setName(int exactUpTo) {
  name = "IAQint: exact up to degree " + std::to_string(exactUpTo);
}

template class QintT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<>
void QtriSymT<double, SpaceDimension, TimeDimension>::setName(int exactUpTo) {
  name = "Qtri: exact up to degree " + std::to_string(exactUpTo);
}

template<typename T, int sDim, int tDim>
QtriSymT<T, sDim, tDim>::QtriSymT(int exactUpTo) {
  QtriSymT_data<T> qData = GetQtriSymT_data<T>(exactUpTo);
  this->resize(qData.Size());
  int cnt = 0;
  if (qData.d0.size() > 0) { initDataPair(cnt++, 1.0 / T(3.0), 1.0 / T(3.0), qData.d0[0].weight); }
  for (int i = 0; i < qData.d1.size(); ++i) {
    initDataPair(cnt++, qData.d1[i].A, qData.d1[i].B, qData.d1[i].weight);
    initDataPair(cnt++, qData.d1[i].B, qData.d1[i].A, qData.d1[i].weight);
    initDataPair(cnt++, qData.d1[i].B, qData.d1[i].B, qData.d1[i].weight);
  }
  for (int i = 0; i < qData.d2.size(); ++i) {
    initDataPair(cnt++, qData.d2[i].A, qData.d2[i].B, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].A, qData.d2[i].C, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].B, qData.d2[i].A, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].B, qData.d2[i].C, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].C, qData.d2[i].A, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].C, qData.d2[i].B, qData.d2[i].weight);
  }
  setName(exactUpTo);
}

template class QtriSymT<double, SpaceDimension, TimeDimension>;

#ifdef BUILD_IA

template<>
void QtriSymT<IAInterval, SpaceDimension, TimeDimension>::setName(int exactUpTo) {
  name = "IAQtri: exact up to degree " + std::to_string(exactUpTo);
}

template class QtriSymT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<>
QquadT<double, SpaceDimension, TimeDimension>::QquadT(int exactUpTo, int exactUpToY) :
    QTensorT<double, SpaceDimension, TimeDimension>(createQuads(exactUpTo, exactUpToY)) {
  this->name = "Qquad: exact up to degree " + std::to_string(exactUpTo) + " in x and "
               + std::to_string(exactUpToY) + " in y";
}

#ifdef BUILD_IA

template<>
QquadT<IAInterval, SpaceDimension, TimeDimension>::QquadT(int exactUpTo, int exactUpToY) :
    QTensorT<IAInterval, SpaceDimension, TimeDimension>(createQuads(exactUpTo, exactUpToY)) {
  this->name = "IAQquad: exact up to degree " + std::to_string(exactUpTo) + " in x and "
               + std::to_string(exactUpToY) + " in y";
}

#endif

template<>
void QtetSymT<double, SpaceDimension, TimeDimension>::setName(int exactUpTo) {
  name = "Qtet: exact up to degree " + std::to_string(exactUpTo);
}

template<typename T, int sDim, int tDim>
QtetSymT<T, sDim, tDim>::QtetSymT(int exactUpTo) {
  QtetSymT_data<T> qData = GetQtetSymT_data<T>(exactUpTo);
  this->resize(qData.Size());
  int cnt = 0;
  if (qData.d0.size() > 0) { initDataPair(cnt++, 0.25, 0.25, 0.25, qData.d0[0].weight); }
  for (int i = 0; i < qData.d1.size(); ++i) {
    initDataPair(cnt++, qData.d1[i].A, qData.d1[i].B, qData.d1[i].B, qData.d1[i].weight);
    initDataPair(cnt++, qData.d1[i].B, qData.d1[i].A, qData.d1[i].B, qData.d1[i].weight);
    initDataPair(cnt++, qData.d1[i].B, qData.d1[i].B, qData.d1[i].A, qData.d1[i].weight);
    initDataPair(cnt++, qData.d1[i].B, qData.d1[i].B, qData.d1[i].B, qData.d1[i].weight);
  }
  for (int i = 0; i < qData.d2.size(); ++i) {
    initDataPair(cnt++, qData.d2[i].A, qData.d2[i].A, qData.d2[i].B, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].A, qData.d2[i].B, qData.d2[i].A, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].A, qData.d2[i].B, qData.d2[i].B, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].B, qData.d2[i].A, qData.d2[i].A, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].B, qData.d2[i].A, qData.d2[i].B, qData.d2[i].weight);
    initDataPair(cnt++, qData.d2[i].B, qData.d2[i].B, qData.d2[i].A, qData.d2[i].weight);
  }
  for (int i = 0; i < qData.d3.size(); ++i) {
    initDataPair(cnt++, qData.d3[i].A, qData.d3[i].B, qData.d3[i].B, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].A, qData.d3[i].B, qData.d3[i].C, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].A, qData.d3[i].C, qData.d3[i].B, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].B, qData.d3[i].A, qData.d3[i].B, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].B, qData.d3[i].A, qData.d3[i].C, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].B, qData.d3[i].B, qData.d3[i].A, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].B, qData.d3[i].B, qData.d3[i].C, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].B, qData.d3[i].C, qData.d3[i].A, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].B, qData.d3[i].C, qData.d3[i].B, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].C, qData.d3[i].A, qData.d3[i].B, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].C, qData.d3[i].B, qData.d3[i].A, qData.d3[i].weight);
    initDataPair(cnt++, qData.d3[i].C, qData.d3[i].B, qData.d3[i].B, qData.d3[i].weight);
  }
  for (int i = 0; i < qData.d4.size(); ++i) {
    initDataPair(cnt++, qData.d4[i].A, qData.d4[i].B, qData.d4[i].C, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].A, qData.d4[i].B, qData.d4[i].D, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].A, qData.d4[i].C, qData.d4[i].B, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].A, qData.d4[i].C, qData.d4[i].D, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].A, qData.d4[i].D, qData.d4[i].B, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].A, qData.d4[i].D, qData.d4[i].C, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].B, qData.d4[i].A, qData.d4[i].C, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].B, qData.d4[i].A, qData.d4[i].D, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].B, qData.d4[i].C, qData.d4[i].A, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].B, qData.d4[i].C, qData.d4[i].D, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].B, qData.d4[i].D, qData.d4[i].A, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].B, qData.d4[i].D, qData.d4[i].C, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].C, qData.d4[i].A, qData.d4[i].B, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].C, qData.d4[i].A, qData.d4[i].D, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].C, qData.d4[i].B, qData.d4[i].A, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].C, qData.d4[i].B, qData.d4[i].D, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].C, qData.d4[i].D, qData.d4[i].A, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].C, qData.d4[i].D, qData.d4[i].B, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].D, qData.d4[i].A, qData.d4[i].B, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].D, qData.d4[i].A, qData.d4[i].C, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].D, qData.d4[i].B, qData.d4[i].A, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].D, qData.d4[i].B, qData.d4[i].C, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].D, qData.d4[i].C, qData.d4[i].A, qData.d4[i].weight);
    initDataPair(cnt++, qData.d4[i].D, qData.d4[i].C, qData.d4[i].B, qData.d4[i].weight);
  }
  setName(exactUpTo);
}

template class QtetSymT<double, SpaceDimension, TimeDimension>;

#ifdef BUILD_IA

template<>
void QtetSymT<IAInterval, SpaceDimension, TimeDimension>::setName(int exactUpTo) {
  name = "IAQtet: exact up to degree " + std::to_string(exactUpTo);
}

template class QtetSymT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<>
QpriSymT<double, SpaceDimension, TimeDimension>::QpriSymT(int exactUpToTriangle, int exactUpToZ) :
    QTensorT<double, SpaceDimension, TimeDimension>(createQuads(exactUpToTriangle, exactUpToZ),
                                                    true) {
  this->name = "Qpri: exact up to degree " + std::to_string(exactUpToTriangle) + " on triangle and "
               + std::to_string(exactUpToZ) + " in z";
}

#ifdef BUILD_IA

template<>
QpriSymT<IAInterval, SpaceDimension, TimeDimension>::QpriSymT(int exactUpToTriangle,
                                                              int exactUpToZ) :
    QTensorT<IAInterval, SpaceDimension, TimeDimension>(createQuads(exactUpToTriangle, exactUpToZ),
                                                        true) {
  this->name = "IAQpri: exact up to degree " + std::to_string(exactUpToTriangle)
               + " on triangle and " + std::to_string(exactUpToZ) + " in z";
}

#endif

template<>
QhexT<double, SpaceDimension, TimeDimension>::QhexT(int exactUpTo, int exactUpToY, int exactUpToZ) :
    QTensorT<double, SpaceDimension, TimeDimension>(
        createQuads(exactUpTo, exactUpToY, exactUpToZ)) {
  this->name = "Qhex: exact up to degree " + std::to_string(exactUpTo) + " in x and "
               + std::to_string(exactUpToY) + " in y and " + std::to_string(exactUpToZ) + " in z";
}

#ifdef BUILD_IA

template<>
QhexT<IAInterval, SpaceDimension, TimeDimension>::QhexT(int exactUpTo, int exactUpToY,
                                                        int exactUpToZ) :
    QTensorT<IAInterval, SpaceDimension, TimeDimension>(
        createQuads(exactUpTo, exactUpToY, exactUpToZ)) {
  this->name = "IAQhex: exact up to degree " + std::to_string(exactUpTo) + " in x and "
               + std::to_string(exactUpToY) + " in y and " + std::to_string(exactUpToZ) + " in z";
}

#endif
