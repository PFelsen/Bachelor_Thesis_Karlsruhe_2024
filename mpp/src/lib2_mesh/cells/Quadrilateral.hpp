#ifndef QUADRILATERAL_H
#define QUADRILATERAL_H

#include "ICell.hpp"

// ----------------------------------------------------
//    QUADRILATERAL
// ----------------------------------------------------
class Quadrilateral : public Cell {
  template<typename T>
  PointT<T> globalToLocal(const PointT<T> &z) const {
#ifndef AFFINE_TRANSFORMATION
    // https://math.stackexchange.com/questions/3069262/functions-that-map-a-quadrilateral-to-the-unit-square
    PointT<T> a = corners[0];
    PointT<T> b = -corners[0] + corners[1];
    PointT<T> c = -corners[0] + corners[3];
    PointT<T> d = corners[0] - corners[1] + corners[2] - corners[3];

    if (norm(d) < 1e-14) {
      auto Tran = GetTransformation(z);
      return Tran.ApplyJinv(z - corners[0]);
    }
    auto A = d[0] * c[1] - c[0] * d[1];
    auto B = d[0] * (a[1] - z[1]) + d[1] * (z[0] - a[0]) + b[0] * c[1] - c[0] * b[1];
    auto C = b[0] * (a[1] - z[1]) + b[1] * (z[0] - a[0]);
    PointT<T> res;
    if (abs(A) < 1e-14) {
      if (abs(B) < 1e-14) {
        THROW("Quadrilateral::globalToLocal A == B == 0")
      } else {
        res[1] = -C / B;
      }
    } else {
      auto p = B / A;
      auto q = C / A;
      auto w = p * p / 4 - q;
      if (w < 0.0) { THROW("Quadrilateral::globalToLocal Discriminant negative.") }
      res[1] = -p / 2 + sqrt(w);
      if (res[1] < 0 || 1 < res[1]) { res[1] = -p / 2 - sqrt(w); }
    }
    if (abs(b[0] + d[0] * res[1]) < 1e-14) {
      // THROW("Quadrilateral::globalToLocal dividing by zero.")
      auto Tran = GetTransformation(z);
      return Tran.ApplyJinv(z - corners[0]);
    }
    res[0] = (z[0] - a[0] - c[0] * res[1]) / (b[0] + d[0] * res[1]);
    return res;
#endif
    auto Tran = GetTransformation(z);
    return Tran.ApplyJinv(z - corners[0]);
  }

  template<typename T>
  PointT<T> localToGlobal(const PointT<T> &local) const {
    PointT<T> v = (1.0 - local[0]) * corners[0] + local[0] * corners[1];
    PointT<T> w = (1.0 - local[0]) * corners[3] + local[0] * corners[2];
    return (1.0 - local[1]) * v + local[1] * w;
  }

  template<typename T>
  PointT<T> faceLocalToGlobal(int face, const PointT<T> &localOnFace) const {
    return (1 - localOnFace[0]) * FaceCorner(face, 0) + localOnFace[0] * FaceCorner(face, 1);
  }

  template<typename T>
  PointT<T> faceLocalToLocal(int face, const PointT<T> &localOnFace) const {
    return (1 - localOnFace[0]) * LocalFaceCorner(face, 0)
           + localOnFace[0] * LocalFaceCorner(face, 1);
  }

  template<typename T>
  PointT<T> faceLocalToLocalOriented(int faceId, const PointT<T> &localOnFace) const {
    if (FaceCorner(faceId, 0) < FaceCorner(faceId, 1))
      return (1 - localOnFace[0]) * LocalFaceCorner(faceId, 0)
             + localOnFace[0] * LocalFaceCorner(faceId, 1);
    return (1 - localOnFace[0]) * LocalFaceCorner(faceId, 1)
           + localOnFace[0] * LocalFaceCorner(faceId, 0);
  }

  template<typename T, int sDim, int tDim>
  TransformationT<T, sDim> getTransformation(const PointT<T, sDim, tDim> &z) const {
    vector<PointT<T, sDim, tDim>> y(2);
    for (int k = 0; k < y.size(); ++k) {
      y[k] = PointT<T, sDim, tDim>(-1 + z[1], -1 + z[0]) * corners[0][k]
             + PointT<T, sDim, tDim>(1 - z[1], -z[0]) * corners[1][k]
             + PointT<T, sDim, tDim>(z[1], z[0]) * corners[2][k]
             + PointT<T, sDim, tDim>(-z[1], 1 - z[0]) * corners[3][k];
    }
    return TransformationT<T, sDim>(y);
  }
public:
  Quadrilateral(const vector<Point> &z, short sd);

  static const Quadrilateral *reference;

  static const Point *LocalCornersQuad();

  static const Point *LocalEdgesQuad();

  static const Point *LocalFaceNormalsQuad();

  static const Point *LocalCenterQuad();

  static const vector<Rule> *Rquad;

  static const vector<Rule> *Rquad_barycentric;

  virtual CELLTYPE Type() const override;

  CELLTYPE ReferenceType() const override;

  int Corners() const override;

  const Point &LocalCorner(int i) const override;

  const Point &LocalEdge(int i) const override;

  const Point &LocalFace(int i) const override;

  const Point &LocalFaceNormal(int i) const override;

  const Point &LocalCenter() const override;

  virtual Point FaceLocalToGlobal(int face, const Point &) const override;

  virtual Point FaceLocalToLocal(int faceId, const Point &localOnFace) const override;

  virtual Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const override;

  virtual Point LocalToGlobal(const Point &z) const override;

  virtual Point GlobalToLocal(const Point &point) const override;

#ifdef BUILD_IA
  IAPoint GlobalToLocal(const IAPoint &point) const override;

  IAPoint LocalToGlobal(const IAPoint &z) const override;

  IAPoint FaceLocalToGlobal(int face, const IAPoint &) const override;

  IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const override;

  IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const override;

  IATransformation GetTransformation(const IAPoint &x) const override;

  IAInterval LocalFaceAreaIA(int faceId) const override;

  static const IAPoint *LocalFaceNormalsIAQuad();

  const IAPoint &LocalFaceNormalIA(int faceId) const override;
#endif

  bool PointInCell(const Point &z) const override;

  virtual Transformation GetTransformation(const Point &x) const override;

  int Edges() const override;

  Point Edge(int i) const override;

  short edgecorner(unsigned short i, unsigned short j) const override;

  int Faces() const override;

  Point Face(int i) const override;

  CELLTYPE FaceType(int face) const override;

  short FaceCorners(int i) const override;

  short facecorner(unsigned short i, unsigned short j) const override;

  short FaceEdges(int i) const override;

  short faceedge(int i, int j) const override;

  double FaceArea(int i) const override;

  double LocalFaceArea(int i) const override;

  void GetRefinePoints(std::vector<Point> &z) const override;

  virtual const vector<Rule> &Refine(vector<Point> &a) const override;

  void GetRefineBarycentricPoints(std::vector<Point> &z, const Point &shiftCenter) const override;

  const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                             const Point &shiftCenter) const override;

  int dim() const override;

  bool plane() const override;

  static const Quadrilateral *ReferenceQuadrilateral();

  const Cell *ReferenceCell() const override;

  int Children() const override;

  Point Child(int i) const override;
};

class Quadrilateral2 : public Quadrilateral {
public:
  static const vector<Rule> *R2quad;

  Quadrilateral2(const vector<Point> &z, short sd);

  Point FaceLocalToGlobal(int face, const Point &) const override{
      THROW("Check if FaceLocalToGlobal needs to be implemented")}

  Point FaceLocalToLocal(int face, const Point &) const override{
      THROW("Check if FaceLocalToLocal needs to be implemented")}

  Point FaceLocalToLocalOriented(int face, const Point &) const override{
      THROW("Check if FaceLocalToLocalOriented needs to be implemented")}

  Point LocalToGlobal(const Point &z) const override;

  Transformation GetTransformation(const Point &x) const override;

  CELLTYPE Type() const override;

  void GetRefinePoints(std::vector<Point> &a) const override;

  const vector<Rule> &Refine(vector<Point> &z) const override;

  const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                             const Point &shiftCenter) const override {
    THROW("Barycentric refinement not implemented for cell type " + std::to_string(Type()))
  }

  static const Quadrilateral2 *ReferenceQuadrilateral2();

#ifdef BUILD_IA
  IAPoint GlobalToLocal(const IAPoint &point) const override{
      THROW("GlobalToLocal not implemented for IA types")}

  IAPoint LocalToGlobal(const IAPoint &z) const override{
      THROW("LocalToGlobal not implemented for IA types")}

  IAPoint FaceLocalToGlobal(int face, const IAPoint &) const override{
      THROW("FaceLocalToGlobal not implemented for IA types")}

  IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const override{
      THROW("FaceLocalToLocal not implemented for IA types")}

  IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const override{
      THROW("FaceLocalToLocalOriented not implemented for IA types")}

  IATransformation GetTransformation(const IAPoint &x) const override {
    THROW("GetTransformation not implemented for IA types")
  }
#endif
};

#endif // QUADRILATERAL_H
