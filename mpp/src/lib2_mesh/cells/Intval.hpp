#ifndef INTVAL_H
#define INTVAL_H

#include "ICell.hpp"

// ----------------------------------------------------
//    INTERVAL
// ----------------------------------------------------

class Intval : public Cell {
  template<typename T>
  PointT<T> globalToLocal(const PointT<T> &z) const {
    PointT<T> c1(corners[1]); // important for interval arithmetic!
    return PointT<T>((z[0] - corners[0][0]) / (c1[0] - corners[0][0]), 0.0, 0.0);
  }

  template<typename T>
  PointT<T> localToGlobal(const PointT<T> &local) const {
    return (1 - local[0]) * corners[0] + local[0] * corners[1];
  }

  template<typename T>
  PointT<T> faceLocalToGlobal(int face, const PointT<T> &localOnFace) const {
    return localOnFace[0] * FaceCorner(face, 0);
  }

  template<typename T>
  PointT<T> faceLocalToLocal(int face, const PointT<T> &localOnFace) const {
    return localOnFace[0] * LocalFaceCorner(face, 0);
  }

  template<typename T>
  PointT<T> faceLocalToLocalOriented(int face, const PointT<T> &localOnFace) const {
    return localOnFace[0] * LocalFaceCorner(face, 0);
  }

  template<typename T>
  TransformationT<T> getTransformation(const PointT<T> &z) const {
    vector<PointT<T>> y(1);
    y[0] = PointT<T>(-1.0, 0.0) * corners[0][0] + PointT<T>(1.0, 0.0) * corners[1][0];
    return TransformationT<T>(y);
  }
public:
  static const Point *LocalCornersInt();

  static const Point *LocalEdgesInt();

  static const Point *LocalFaceNormalsInt();

  static const Point *LocalCenterInt();

  static const vector<Rule> *Rint;

  Intval(const vector<Point> &x, int sd);

  CELLTYPE Type() const override;

  CELLTYPE ReferenceType() const override;

  int Corners() const override;

  const Point &LocalCorner(int i) const override;

  const Point &LocalEdge(int i) const override;

  const Point &LocalFace(int i) const override;

  const Point &LocalFaceNormal(int i) const override;

  const Point &LocalCenter() const override;

  Point FaceLocalToGlobal(int face, const Point &) const override;

  Point FaceLocalToLocal(int faceId, const Point &localOnFace) const override;

  Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const override;

  Point GlobalToLocal(const Point &point) const override;

  Point LocalToGlobal(const Point &z) const override;

  Transformation GetTransformation(const Point &x) const override;

#ifdef BUILD_IA
  IAPoint GlobalToLocal(const IAPoint &point) const override;

  IAPoint LocalToGlobal(const IAPoint &z) const override;

  IAPoint FaceLocalToGlobal(int face, const IAPoint &) const override;

  IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const override;

  IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const override;

  IATransformation GetTransformation(const IAPoint &x) const override;

  IAInterval LocalFaceAreaIA(int faceId) const override;

  static const IAPoint *LocalFaceNormalsIAInt();

  const IAPoint &LocalFaceNormalIA(int faceId) const override;
#endif

  int Edges() const override;

  int EdgeCorners(int edgeId) const override { return 1; }

  Point Edge(int i) const override;

  short edgecorner(unsigned short edge, unsigned short corner) const override;

  int Faces() const override;

  Point Face(int i) const override;

  CELLTYPE FaceType(int face) const override;

  short FaceCorners(int i) const override;

  short facecorner(unsigned short face, unsigned short corner) const override;

  short FaceEdges(int i) const override;

  short faceedge(int i, int j) const override;

  double FaceArea(int i) const override;

  double LocalFaceArea(int i) const override;

  void GetRefinePoints(std::vector<Point> &a) const override;

  const vector<Rule> &Refine(vector<Point> &a) const override;

  int dim() const override;

  bool plane() const override;

  static const Intval *ReferenceIntval();

  const Cell *ReferenceCell() const override;

  int Children() const override;

  Point Child(int i) const override;

  virtual bool PointInCell(const Point &z) const override;
};

#endif // INTVAL_H
