#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include <utility>
#include "ICell.hpp"

// ----------------------------------------------------
//    TETRAHEDRON
// ----------------------------------------------------
class Tetrahedron : public Cell {
  template<typename POINT>
  POINT globalToLocal(const POINT &z) const {
    auto T = GetTransformation(z);
    return T.ApplyJinv(z - corners[0]);
  }

  template<typename POINT>
  POINT localToGlobal(const POINT &local) const {
    return (1 - local[0] - local[1] - local[2]) * corners[0] + local[0] * corners[1]
           + local[1] * corners[2] + local[2] * corners[3];
  }

  template<typename POINT>
  POINT faceLocalToGlobal(int face, const POINT &localOnFace) const {
    return (1 - localOnFace[0] - localOnFace[1]) * FaceCorner(face, 0)
           + localOnFace[0] * FaceCorner(face, 1) + localOnFace[1] * FaceCorner(face, 2);
  }

  template<typename POINT>
  POINT faceLocalToLocal(int faceId, const POINT &localOnFace) const {
    return (1 - localOnFace[0] - localOnFace[1]) * LocalFaceCorner(faceId, 0)
           + localOnFace[0] * LocalFaceCorner(faceId, 1)
           + localOnFace[1] * LocalFaceCorner(faceId, 2);
  }

  template<typename POINT>
  POINT faceLocalToLocalOriented(int faceId, const POINT &localOnFace) const {
    const Point &corner0 = Corner(facecorner(faceId, 0));
    const Point &corner1 = Corner(facecorner(faceId, 1));
    const Point &corner2 = Corner(facecorner(faceId, 2));
    const Point &localCorner0 = LocalCorner(facecorner(faceId, 0));
    const Point &localCorner1 = LocalCorner(facecorner(faceId, 1));
    const Point &localCorner2 = LocalCorner(facecorner(faceId, 2));
    if ((corner0 < corner1) && (corner1 < corner2))
      return (1 - localOnFace[0] - localOnFace[1]) * localCorner0 + localOnFace[0] * localCorner1
             + localOnFace[1] * localCorner2;
    else if ((corner0 < corner2) && (corner2 < corner1))
      return (1 - localOnFace[0] - localOnFace[1]) * localCorner0 + localOnFace[0] * localCorner2
             + localOnFace[1] * localCorner1;
    else if ((corner1 < corner2) && (corner2 < corner0))
      return (1 - localOnFace[0] - localOnFace[1]) * localCorner1 + localOnFace[0] * localCorner2
             + localOnFace[1] * localCorner0;
    else if ((corner1 < corner0) && (corner0 < corner2))
      return (1 - localOnFace[0] - localOnFace[1]) * localCorner1 + localOnFace[0] * localCorner0
             + localOnFace[1] * localCorner2;
    else if ((corner2 < corner0) && (corner0 < corner1))
      return (1 - localOnFace[0] - localOnFace[1]) * localCorner2 + localOnFace[0] * localCorner0
             + localOnFace[1] * localCorner1;
    else if ((corner2 < corner1) && (corner1 < corner0))
      return (1 - localOnFace[0] - localOnFace[1]) * localCorner2 + localOnFace[0] * localCorner1
             + localOnFace[1] * localCorner0;
    THROW("Error in FaceLocalToLocal")
  }

  template<typename T, int sDim, int tDim>
  TransformationT<T, sDim> getTransformation(const PointT<T, sDim, tDim> &z) const {
    vector<PointT<T, sDim, tDim>> y(3);
    for (int k = 0; k < y.size(); ++k) {
      y[k] = PointT<T, sDim, tDim>(-1.0, -1.0, -1.0) * corners[0][k]
             + PointT<T, sDim, tDim>(1.0, 0.0, 0.0) * corners[1][k]
             + PointT<T, sDim, tDim>(0.0, 1.0, 0.0) * corners[2][k]
             + PointT<T, sDim, tDim>(0.0, 0.0, 1.0) * corners[3][k];
    }
    return TransformationT<T, sDim>(y);
  }
public:
  static const Point *LocalCornersTet();

  static const Point *LocalEdgesTet();

  static const Point *LocalFacesTet();

  static const Point *LocalFaceNormalsTet();

  static const Point *LocalCenterTet();

  static const vector<Rule> *Rtet;
  static const vector<Rule> *Rtet_1;
  static const vector<Rule> *Rtet_2;

  static const vector<Rule> *Rtet_barycentric;

  Tetrahedron(const vector<Point> &z, int sd);

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

  Point LocalToGlobal(const Point &z) const override;

  Point GlobalToLocal(const Point &point) const override;

  Transformation GetTransformation(const Point &x) const override;

#ifdef BUILD_IA
  IAPoint GlobalToLocal(const IAPoint &point) const override;

  IAPoint LocalToGlobal(const IAPoint &z) const override;

  IAPoint FaceLocalToGlobal(int face, const IAPoint &) const override;

  IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const override;

  IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const override;

  IATransformation GetTransformation(const IAPoint &x) const override;

  IAInterval LocalFaceAreaIA(int faceId) const override;

  static const IAPoint *LocalFaceNormalsIATet();

  const IAPoint &LocalFaceNormalIA(int faceId) const override;
#endif

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

  double LocalFaceArea(int i) const override;

  double GetVolume() const override;

  VectorField GetNablaNi(int i) const override;

  double FaceArea(int i) const override;

  void GetRefinePoints(std::vector<Point> &z) const override;

  const vector<Rule> &Refine(vector<Point> &z) const override;

  void GetRefineBarycentricPoints(std::vector<Point> &z, const Point &shiftCenter) const override;

  const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                             const Point &shiftCenter) const override;

  int dim() const override;

  bool plane() const override;

  static const Tetrahedron *ReferenceTetrahedron();

  const Cell *ReferenceCell() const override;

  int Children() const override;

  Point Child(int i) const override;
};

#endif // TETRAHEDRON_H
