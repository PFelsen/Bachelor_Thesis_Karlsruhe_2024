#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <string>
#include "Cell.hpp"
#include "Tensor.hpp"

template<typename T>
class ShapeValues : public vector<vector<T>> {
public:
  ShapeValues() = default;

  ShapeValues(int quadraturePoints, int nodalPoints) :
      vector<vector<T>>(quadraturePoints, vector<T>(nodalPoints)) {}

  void resize(int quadraturePoints, int nodalPoints) {
    std::vector<vector<T>>::resize(quadraturePoints);
    for (int q = 0; q < quadraturePoints; ++q)
      (*this)[q].resize(nodalPoints);
  }
};

template<typename T>
class FaceShapeValues : public vector<vector<vector<T>>> {
public:
  using vector<vector<vector<T>>>::resize;

  FaceShapeValues() = default;

  FaceShapeValues(int faces, int quadraturePoints, int nodalPoints) :
      vector<vector<vector<T>>>(faces,
                                vector<vector<T>>(quadraturePoints, vector<T>(nodalPoints))) {}

  void resize(int faces, int quadraturePoints, int nodalPoints) {
    std::vector<vector<vector<T>>>::resize(faces);
    for (int f = 0; f < faces; ++f) {
      (*this)[f].resize(quadraturePoints);
      for (int q = 0; q < quadraturePoints; ++q)
        (*this)[f][q].resize(nodalPoints);
    }
  }
};

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ShapeT {
protected:
  ShapeValues<T> localValues = {};
  ShapeValues<VectorFieldT<T, sDim>> localGradient = {};
  ShapeValues<SymTensorT<T, sDim>> localHessian = {};
  ShapeValues<VectorFieldT<T, sDim>> localVector = {};
  ShapeValues<VectorFieldT<T, sDim>> localCurl = {};
  ShapeValues<T> localDivs = {};
  ShapeValues<T> localLaplaces = {};

  FaceShapeValues<T> localFaceValues;

  int numNodalPoints;
public:
  ShapeT(int numNodalPoints) : numNodalPoints(numNodalPoints) {}

  virtual ~ShapeT() {}

  virtual const std::string Name() const = 0;

  int size() const { return numNodalPoints; }

  const ShapeValues<T> &values() const { return localValues; }

  T operator()(int q, int i) const { return localValues[q][i]; }

  virtual T operator()(const PointT<T, sDim, tDim> &, int) const {THROW("Not implemented!")}

  T
  operator()(int fid, int q, int i) const {
    return localFaceValues[fid][q][i];
  }

  VectorFieldT<T, sDim> LocalGradient(int q, int i) const { return localGradient[q][i]; }

  virtual VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &,
                                              int) const {THROW("Not implemented!")}

  SymTensorT<T, sDim> LocalHessian(int q, int i) const {
    return localHessian[q][i];
  }

  virtual SymTensorT<T, sDim> LocalHessian(const PointT<T, sDim, tDim> &,
                                           int) const {THROW("Not implemented!")}

  VectorFieldT<T, sDim> LocalVector(int q, int i) const {
    return localVector[q][i];
  }

  virtual VectorFieldT<T, sDim> LocalVector(const PointT<T, sDim, tDim> &,
                                            int) const {THROW("Not implemented!")}

  VectorFieldT<T, sDim> LocalCurl(int q, int i) const {
    return localCurl[q][i];
  }

  virtual VectorFieldT<T, sDim> LocalCurl(const PointT<T, sDim, tDim> &,
                                          int) const {THROW("Not implemented!")}

  T LocalDiv(int q, int i) const {
    return localDivs[q][i];
  }

  virtual T LocalDiv(const PointT<T, sDim, tDim> &, int) const {THROW("Not implemented!")}

  T LocalLaplace(int q, int i) const {
    return localLaplaces[q][i];
  }

  virtual T LocalLaplace(const PointT<T, sDim, tDim> &, int) const { THROW("Not implemented!") }

  virtual void fillValues(const vector<PointT<T, sDim, tDim>> &localQuadraturePoints) = 0;

  virtual void
  fillFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadraturePoints) = 0;

  /**
   * Sets (in contrast to DoF::GetNodalPoints) the true nodal points on a cell
   */
  virtual void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &points) const {
    THROW("Not implemented!")
  }

  friend std::ostream &operator<<(std::ostream &s, const ShapeT<T, sDim, tDim> &S) {
    return s << S.Name();
  }
};

typedef ShapeT<> Shape;

#endif // of #ifndef _SHAPE_H_
