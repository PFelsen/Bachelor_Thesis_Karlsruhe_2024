#ifndef ELEMENT_H
#define ELEMENT_H

#include "ElementData.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class BasicElementT : public rows {
protected:
  const VectorMatrixBase &base;
  const Cell &c; // TODO: move cell to elementData?
  std::shared_ptr<ElementDataT<TT, sDim, tDim>> elementData = nullptr;

  const row &r(int i) const { return (*this)[i]; }
public:
  BasicElementT(const VectorMatrixBase &base, const Cell &c) :
      rows(base.GetMatrixGraph(), c), base(base), c(c) {
    elementData =
        std::make_shared<ElementDataT<TT, sDim, tDim>>(c,
                                                       base.GetDiscT<TT, sDim, tDim>().GetQuad(
                                                           this->c),
                                                       base.GetMatrixGraph().GetMesh());
  }

  BasicElementT(const VectorMatrixBase &base, const Cell &c, int face) :
      rows(base.GetMatrixGraph(), c), base(base), c(c) {
    elementData =
        std::make_shared<FaceElementDataT<TT, sDim, tDim>>(c, face,
                                                           base.GetDiscT<TT, sDim, tDim>()
                                                               .GetFaceQuad(c),
                                                           base.GetMatrixGraph().GetMesh());
  }

  BasicElementT(const VectorMatrixBase &base, const Cell &c, int face, int n) :
      rows(n, base.GetMatrixGraph(), c), base(base), c(c) {
    elementData =
        std::make_shared<FaceElementDataT<TT, sDim, tDim>>(c, face,
                                                           base.GetDiscT<TT, sDim, tDim>()
                                                               .GetFaceQuad(c),
                                                           base.GetMatrixGraph().GetMesh());
  }

  BasicElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
      base(base), c(baseElement.c), rows(base.GetMatrixGraph(), baseElement.c),
      elementData(baseElement.elementData) {}

  rows GetRowsOnFace(int face) const { return rows(base.GetMatrixGraph(), c, face); }

  /**
   * Function to iterate over quadrature data id and weight:
   * Usage:
   * 1. for (auto q : element.Quad()) { -> q.id and q.weight available }
   * 2. for (auto [id, weight] : element.Quad()) {}
   */
  typename QuadPairT<TT, sDim, tDim>::Iter Quad() const {
    return typename QuadPairT<TT, sDim, tDim>::Iter{*elementData};
  }

  /**
   * Function to iterate over quadrature data id, weight and quad point:
   * Usage:
   * 1. for (auto q : element.QuadWithPoint()) { -> q.id, q.weight and q.point available }
   * 2. for (auto [id, weight, point] : element.QuadWithPoint()) {}
   */
  typename QuadTripleT<TT, sDim, tDim>::Iter QuadWithPoint() const {
    return typename QuadTripleT<TT, sDim, tDim>::Iter{*elementData};
  }

  /**
   * Function to iterate over quadrature data id, weight, quad point and local quad point:
   * Usage:
   * 1. for (auto q : element.QuadWithPointAndLocalPoint()) { -> q.id, q.weight and q.point,
   * q.localPoint available }
   * 2. for (auto [id, weight, point, localPoint] : element.QuadWithPointAndLocalPoint()) {}
   */
  typename QuadQuadrupleT<TT, sDim, tDim>::Iter QuadWithPointAndLocalPoint() const {
    return typename QuadQuadrupleT<TT, sDim, tDim>::Iter{*elementData};
  }

  bool Bnd() const { return elementData->bfParts.onBnd(); }

  int Bnd(int i) const { return elementData->bfParts[i]; }

  const VectorMatrixBase &GetVectorMatrixBase() const { return base; }

  const Cell &GetCell() const { return c; }

  int nQ() const { return elementData->quad.size(); }

  TT QWeight(int q) const { return elementData->quadWeights[q]; }

  virtual int NodalPoints() const { return this->size(); }

  virtual const PointT<TT, sDim, tDim> &NodalPoint(int i) const {
    if constexpr (std::is_same<TT, double>::value) return (*this)[i]();
    THROW("Not implemented")
  }

  int NodalPointIndex(const Point &z) const {
    for (int i = 0; i < this->size(); ++i)
      if ((*this)[i]() == z) return i;
    THROW("NodalPoint not found")
  }

  /// Returns the q-th quadrature point in local coordinates of the face
  const PointT<TT, sDim, tDim> &LocalQPoint(int q) const { return elementData->quad.QPoint(q); }

  /// Returns the q-th quadrature point in global coordinates
  const PointT<TT, sDim, tDim> &QPoint(int q) const { return elementData->quadPoints[q]; }

  const TransformationT<TT, sDim> &GetTransformation(int q) const {
    return elementData->GetTransformation(q);
  }

  TransformationT<TT, sDim> GetTransformation(const PointT<TT, sDim, tDim> &local) const {
#ifdef AFFINE_TRANSFORMATION
    return elementData->trafo;
#else
    return c.GetTransformation(local);
#endif
  }

  TransformationT<TT, sDim> GetTransformationCell(const PointT<TT, sDim, tDim> &z) {
    return c.GetTransformation(z);
  }

  PointT<TT, sDim, tDim> ApplyJacobian(int q, const PointT<TT, sDim, tDim> &local) const {
    return this->GetTransformation(q)(local);
  }

  TT GetTransformationDet(int q) const { return this->GetTransformation(q).Det(); }

  /// Returns the q-th quadrature point in local coordinates of the cell
  const PointT<TT, sDim, tDim> &QLocal(int q) const { return elementData->QLocal(q); }

  const PointT<TT, sDim, tDim> &QNormal(int q) const { return elementData->QNormal(q); }

  short Sign() const { return elementData->Sign(); }

  int FaceID() const { return elementData->FaceID(); }

  TT Area() const { return elementData->Area(); }

  virtual ~BasicElementT() = default;
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ElementT : public BasicElementT<TT, sDim, tDim> {
protected:
  const ShapeT<TT, sDim, tDim> &shape;
public:
  ElementT(const VectorMatrixBase &base, const Cell &c) :
      BasicElementT<TT, sDim, tDim>(base, c),
      shape(base.GetDiscT<TT, sDim, tDim>().GetShape(this->c)) {}

  ElementT(const VectorMatrixBase &base, const Cell &c, int face) :
      BasicElementT<TT, sDim, tDim>(base, c, face),
      shape(base.GetDiscT<TT, sDim, tDim>().GetShape(this->c)) {}

  ElementT(const VectorMatrixBase &base, const Cell &c, int face, int n) :
      BasicElementT<TT, sDim, tDim>(base, c, face, n),
      shape(base.GetDiscT<TT, sDim, tDim>().GetShape(this->c, n)) {}

  ElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
      BasicElementT<TT, sDim, tDim>(base, baseElement),
      shape(base.GetDiscT<TT, sDim, tDim>().GetShape(this->c)) {}
};

using Element = ElementT<>;

#include "RTLagrangeDiscretization.hpp"

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedElementT : public BasicElementT<TT, sDim, tDim> {
protected:
  const MixedDoF &mixedDoF;
  const vector<vector<SPData>> &nodalPointData;
  vector<const ShapeT<TT, sDim, tDim> *> shapes;
public:
  MixedElementT(const VectorMatrixBase &base, const Cell &c) :
      BasicElementT<TT, sDim, tDim>(base, c), mixedDoF(MixedDoF::Cast(base.GetShapeDoF())),
      nodalPointData(mixedDoF.StoragePointData(c)), shapes(mixedDoF.NumberOfDoFs()) {
    for (int n = 0; n < shapes.size(); ++n)
      shapes[n] = &base.GetDiscT<TT, sDim, tDim>().GetShape(this->c, n);
  }

  MixedElementT(const VectorMatrixBase &base, const Cell &c, int face) :
      BasicElementT<TT, sDim, tDim>(base, c, face), mixedDoF(MixedDoF::Cast(base.GetShapeDoF())),
      nodalPointData(mixedDoF.StoragePointData(c)), shapes(mixedDoF.NumberOfDoFs()) {

    for (int n = 0; n < shapes.size(); ++n)
      shapes[n] = &base.GetDiscT<TT, sDim, tDim>().GetShape(this->c, n);
  }

  MixedElementT(const VectorMatrixBase &base, const Cell &c, int face, int n) :
      BasicElementT<TT, sDim, tDim>(base, c, face, n), mixedDoF(MixedDoF::Cast(base.GetShapeDoF())),
      nodalPointData(mixedDoF.StoragePointData(c)), shapes(mixedDoF.NumberOfDoFs()) {

    for (int n = 0; n < shapes.size(); ++n)
      shapes[n] = &base.GetDiscT<TT, sDim, tDim>().GetShape(this->c, n);
  }

  MixedElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement) :
      BasicElementT<TT, sDim, tDim>(base, baseElement),
      mixedDoF(MixedDoF::Cast(base.GetShapeDoF())),
      nodalPointData(mixedDoF.StoragePointData(this->c)), shapes(mixedDoF.NumberOfDoFs()) {
    for (int n = 0; n < shapes.size(); ++n)
      shapes[n] = &base.GetDiscT<TT, sDim, tDim>().GetShape(this->c, n);
  }

  rows GetRowsOnFace(int face, int n) const {
    return rows(this->base.GetMatrixGraph(), this->c, face, n);
  }

  const MixedDoF &GetDoF() const { return mixedDoF; }

  int size(int n) const { return shapes[n]->size(); }

  const vector<vector<SPData>> &NodalPointData() const { return nodalPointData; }
};

using MixedElement = MixedElementT<>;

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class IFaceElementT {
  ShapeValues<TT> faceValues{};
protected:
  void ResizeFaceValues(int quadPoints, int nodalPoints) {
    faceValues.resize(quadPoints, nodalPoints);
  }

  void SetFaceValue(int q, int i, TT value) { faceValues[q][i] = value; }
public:
  IFaceElementT(){};

  const ShapeValues<TT> &values() const { return faceValues; }
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class IMixedFaceElementT {
protected:
  vector<ShapeValues<TT>> faceValues{};
protected:
  void SetFaceValue(int q, int n, int i, TT value) { faceValues[n][q][i] = value; }
public:
  IMixedFaceElementT(int dofCount) : faceValues(dofCount){};

  const vector<ShapeValues<TT>> &values() const { return faceValues; }

  const ShapeValues<TT> &values(int n) const { return faceValues[n]; }
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class BasicFaceElementT : public rows {
protected:
  std::shared_ptr<FaceElementDataT<TT, sDim, tDim>> elementData = nullptr;
  const Cell &c;

  const row &r(int i) const { return (*this)[i]; }
public:
  BasicFaceElementT(const VectorMatrixBase &base, const Cell &c, int face) :
      rows(base.GetMatrixGraph(), c, face), c(c) {
    elementData =
        std::make_shared<FaceElementDataT<TT, sDim, tDim>>(c, face,
                                                           base.GetDiscT<TT, sDim, tDim>()
                                                               .GetFaceQuad(c),
                                                           base.GetMatrixGraph().GetMesh());
  }

  BasicFaceElementT(const VectorMatrixBase &base, const Cell &c, int face, int n) :
      rows(base.GetMatrixGraph(), c, face, n), c(c) {
    elementData =
        std::make_shared<FaceElementDataT<TT, sDim, tDim>>(c, face,
                                                           base.GetDiscT<TT, sDim, tDim>()
                                                               .GetFaceQuad(c),
                                                           base.GetMatrixGraph().GetMesh());
  }

  virtual TT Area() const { return elementData->area; }

  const Cell &GetCell() const { return c; }

  int nQ() const { return elementData->quad.size(); }

  virtual int NodalPoints() const { return this->size(); }

  virtual const PointT<TT, sDim, tDim> &NodalPoint(int i) const {
    if constexpr (std::is_same<TT, double>::value) return (*this)[i]();
    THROW("Not implemented")
  }

  virtual ~BasicFaceElementT() = default;

  /// Returns the q-th quadrature point in local coordinates of the face
  const PointT<TT, sDim, tDim> &LocalQPoint(int q) const { return elementData->quad.QPoint(q); }

  /// Returns the q-th quadrature point in local coordinates of the cell
  const PointT<TT, sDim, tDim> &QLocal(int q) const { return elementData->quadPointCellLocals[q]; }

  /// Returns the q-th quadrature point in global coordinates
  const PointT<TT, sDim, tDim> &QPoint(int q) const { return elementData->quadPoints[q]; }

  const PointT<TT, sDim, tDim> &QNormal(int q) const { return elementData->quadNormals[q]; }

  TT QWeight(int q) const { return elementData->quadWeights[q]; }

  const TransformationT<TT, sDim> &GetTransformation(int q) const {
    return elementData->GetTransformation(q);
  }

  TT GetTransformationDet(int q) const { return elementData->GetTransformation(q).Det(); }

  TransformationT<TT, sDim> GetTransformationCell(const PointT<TT, sDim, tDim> &z) {
    return c.GetTransformation(z);
  }

  short Sign() const { return elementData->sign; }

  int FaceID() const { return elementData->faceID; }
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class FaceElementT : public BasicFaceElementT<TT, sDim, tDim> {
protected:
  const ShapeT<TT, sDim, tDim> &shape;
public:
  FaceElementT(const VectorMatrixBase &base, const Cell &c, int face) :
      BasicFaceElementT<TT, sDim, tDim>(base, c, face),
      shape(base.GetDiscT<TT, sDim, tDim>().GetShape(this->c)) {}

  FaceElementT(const VectorMatrixBase &base, const Cell &c, int face, int n) :
      BasicFaceElementT<TT, sDim, tDim>(base, c, face, n),
      shape(base.GetDiscT<TT, sDim, tDim>().GetShape(this->c, n)) {}
};

typedef FaceElementT<> FaceElement;

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class MixedFaceElementT : public BasicFaceElementT<TT, sDim, tDim> {
protected:
  vector<const ShapeT<TT, sDim, tDim> *> shapes;
public:
  MixedFaceElementT(const VectorMatrixBase &base, const Cell &c, int face) :
      BasicFaceElementT<TT, sDim, tDim>(base, c, face), shapes(base.GetDoF().NumberOfDoFs()) {
    for (int m = 0; m < shapes.size(); ++m)
      shapes[m] = &base.GetDiscT<TT, sDim, tDim>().GetShape(c, m);
  }

  MixedFaceElementT(const VectorMatrixBase &base, const Cell &c, int face, int n) :
      BasicFaceElementT<TT, sDim, tDim>(base, c, face, n), shapes(base.GetDoF().NumberOfDoFs()) {
    for (int m = 0; m < shapes.size(); ++m)
      shapes[m] = &base.GetDiscT<TT, sDim, tDim>().GetShape(c, m);
  }

  int sizeShape(int n) const { return shapes[n]->size(); }

  int size(int n) const { THROW("TODO"); }
};


#endif // of #ifndef ELEMENT_H
