#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "ICell.hpp"

// ----------------------------------------------------
//    TRIANGLE
// ----------------------------------------------------
class Triangle : public Cell {
  template<typename T>
  PointT<T> globalToLocal(const PointT<T> &z) const {
    auto tranformation = GetTransformation(z);
    return tranformation.ApplyJinv(z - corners[0]);
  }

  template<typename T>
  PointT<T> localToGlobal(const PointT<T> &local) const {
    return (1 - local[0] - local[1]) * corners[0] + local[0] * corners[1] + local[1] * corners[2];
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

  template<typename T>
  TransformationT<T> getTransformation(const PointT<T> &z) const {
    vector<PointT<T>> y(2);
    for (int k = 0; k < y.size(); ++k) {
      y[k] = PointT<T>(-1.0, -1.0) * corners[0][k] + PointT<T>(1.0, 0.0) * corners[1][k]
             + PointT<T>(0.0, 1.0) * corners[2][k];
    }
    return TransformationT<T>(y);
  }
public:
  static const Point *LocalCornersTri();

  static const Point *LocalEdgesTri();

  static const Point *LocalFaceNormalsTri();

  static const Point *LocalCenterTri();

  static const vector<Rule> *Rtri;

  static const vector<Rule> *Rtri_barycentric;

  Triangle(const vector<Point> &z, short sd);

  CELLTYPE Type() const override;

  CELLTYPE ReferenceType() const override;

  inline int Corners() const override { return 3; }

  const Point &LocalCorner(int i) const override;

  const Point &LocalEdge(int i) const override;

  const Point &LocalFace(int i) const override;

  const Point &LocalFaceNormal(int i) const override;

  const Point &LocalCenter() const override;

  Point LocalToGlobal(const Point &z) const override;

  Point GlobalToLocal(const Point &point) const override;

  Point FaceLocalToGlobal(int face, const Point &) const override;

  Point FaceLocalToLocal(int faceId, const Point &localOnFace) const override;

  Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const override;

  Transformation GetTransformation(const Point &x) const override;

#ifdef BUILD_IA
  IAPoint GlobalToLocal(const IAPoint &point) const override;

  IAPoint LocalToGlobal(const IAPoint &z) const override;

  IAPoint FaceLocalToGlobal(int face, const IAPoint &) const override;

  IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const override;

  IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const override;

  IATransformation GetTransformation(const IAPoint &x) const override;

  IAInterval LocalFaceAreaIA(int faceId) const override;

  static const IAPoint *LocalFaceNormalsIATri();

  const IAPoint &LocalFaceNormalIA(int faceId) const override;
#endif

  int Edges() const override;

  Point Edge(int i) const override;

  short edgecorner(unsigned short i, unsigned short j) const override;

  int Faces() const override;

  inline Point Face(int i) const {
    switch (i) {
    case 0:
      return 0.5 * (Corner(0) + Corner(1));
    case 1:
      return 0.5 * (Corner(1) + Corner(2));
    case 2:
      return 0.5 * (Corner(2) + Corner(0));
    }
    THROW("Not Implemented")
  }

  CELLTYPE FaceType(int face) const override;

  short FaceCorners(int i) const override;

  short facecorner(unsigned short i, unsigned short j) const override;

  short FaceEdges(int i) const override;

  short faceedge(int i, int j) const override;

  double FaceArea(int i) const override;

  double LocalFaceArea(int i) const override;

  void GetRefinePoints(std::vector<Point> &z) const override;

  const vector<Rule> &Refine(vector<Point> &a) const override;

  void GetRefineBarycentricPoints(std::vector<Point> &z, const Point &shiftCenter) const override;

  const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                             const Point &shiftCenter = Origin) const override;

  int dim() const override;

  bool plane() const override;

  static const Triangle *ReferenceTriangle();

  const Cell *ReferenceCell() const override;

  int Children() const override;

  Point Child(int i) const override;

  bool PointInCell(const Point &global) const {
    Point local = GlobalToLocal(global);
    return local[0] >= 0 && local[1] >= 0 && local[0] + local[1] <= 1;
  }
};

#endif // TRIANGLE_H
