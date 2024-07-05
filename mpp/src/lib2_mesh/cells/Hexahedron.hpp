#ifndef HEXAHEDRON_H
#define HEXAHEDRON_H

#include "ICell.hpp"

// ----------------------------------------------------
//    HEXAHEDRON
// ----------------------------------------------------
class Hexahedron : public Cell {
  template<typename POINT>
  POINT globalToLocal(const POINT &z) const {
#ifndef AFFINE_TRANSFORMATION
    Warning("globalToLocal only works for affine transformations")
#endif
        auto T = GetTransformation(z);
    return T.ApplyJinv(z - corners[0]);
  }

  template<typename POINT>
  POINT localToGlobal(const POINT &local) const {
    return (1 - local[0]) * (1 - local[1]) * (1 - local[2]) * corners[0]
           + local[0] * (1 - local[1]) * (1 - local[2]) * corners[1]
           + local[0] * local[1] * (1 - local[2]) * corners[2]
           + (1 - local[0]) * local[1] * (1 - local[2]) * corners[3]
           + (1 - local[0]) * (1 - local[1]) * local[2] * corners[4]
           + local[0] * (1 - local[1]) * local[2] * corners[5]
           + local[0] * local[1] * local[2] * corners[6]
           + (1 - local[0]) * local[1] * local[2] * corners[7];
  }

  template<typename POINT>
  POINT faceLocalToGlobal(int face, const POINT &localOnFace) const {
    return (1 - localOnFace[0]) * (1 - localOnFace[1]) * FaceCorner(face, 0)
           + localOnFace[0] * (1 - localOnFace[1]) * FaceCorner(face, 1)
           + localOnFace[0] * localOnFace[1] * FaceCorner(face, 2)
           + (1 - localOnFace[0]) * localOnFace[1] * FaceCorner(face, 3);
  }

  template<typename POINT>
  POINT faceLocalToLocal(int face, const POINT &localOnFace) const {
    return (1 - localOnFace[0]) * (1 - localOnFace[1]) * LocalFaceCorner(face, 0)
           + localOnFace[0] * (1 - localOnFace[1]) * LocalFaceCorner(face, 1)
           + localOnFace[0] * localOnFace[1] * LocalFaceCorner(face, 2)
           + (1 - localOnFace[0]) * localOnFace[1] * LocalFaceCorner(face, 3);
  }

  template<typename POINT>
  POINT faceLocalToLocalOriented(int faceId, const POINT &localOnFace) const {
    const Point &corner0 = Corner(facecorner(faceId, 0));
    const Point &corner1 = Corner(facecorner(faceId, 1));
    const Point &corner2 = Corner(facecorner(faceId, 2));
    const Point &corner3 = Corner(facecorner(faceId, 3));
    const Point &localCorner0 = LocalCorner(facecorner(faceId, 0));
    const Point &localCorner1 = LocalCorner(facecorner(faceId, 1));
    const Point &localCorner2 = LocalCorner(facecorner(faceId, 2));
    const Point &localCorner3 = LocalCorner(facecorner(faceId, 3));
    if ((corner0 < corner1) && (corner0 < corner2) && (corner0 < corner3)) {
      if (corner1 < corner3)
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner0
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner1
               + localOnFace[0] * localOnFace[1] * localCorner2
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner3;
      else
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner0
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner3
               + localOnFace[0] * localOnFace[1] * localCorner2
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner1;
    } else if ((corner1 < corner0) && (corner1 < corner2) && (corner1 < corner3)) {
      if (corner2 < corner0)
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner1
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner2
               + localOnFace[0] * localOnFace[1] * localCorner3
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner0;
      else
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner1
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner0
               + localOnFace[0] * localOnFace[1] * localCorner3
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner2;
    } else if ((corner2 < corner0) && (corner2 < corner1) && (corner2 < corner3)) {
      if (corner3 < corner1)
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner2
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner3
               + localOnFace[0] * localOnFace[1] * localCorner0
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner1;
      else
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner2
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner1
               + localOnFace[0] * localOnFace[1] * localCorner0
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner3;
    } else if ((corner3 < corner0) && (corner3 < corner2) && (corner3 < corner1)) {
      if (corner0 < corner2)
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner3
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner0
               + localOnFace[0] * localOnFace[1] * localCorner1
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner2;
      else
        return (1 - localOnFace[0]) * (1 - localOnFace[1]) * localCorner3
               + localOnFace[0] * (1 - localOnFace[1]) * localCorner2
               + localOnFace[0] * localOnFace[1] * localCorner1
               + (1 - localOnFace[0]) * localOnFace[1] * localCorner0;
    }
    THROW("Error in FaceGlobalToLocal")
  }

  template<typename T, int sDim, int tDim>
  TransformationT<T, sDim> getTransformation(const PointT<T, sDim, tDim> &z) const {
    vector<PointT<T, sDim, tDim>> y(3);
    for (int k = 0; k < y.size(); ++k) {
      y[k] =
          PointT<T, sDim, tDim>(-(1 - z[1]) * (1 - z[2]), -(1 - z[0]) * (1 - z[2]),
                                -(1 - z[0]) * (1 - z[1]))
              * corners[0][k]
          + PointT<T, sDim, tDim>((1 - z[1]) * (1 - z[2]), -z[0] * (1 - z[2]), -z[0] * (1 - z[1]))
                * corners[1][k]
          + PointT<T, sDim, tDim>(z[1] * (1 - z[2]), z[0] * (1 - z[2]), -z[0] * z[1])
                * corners[2][k]
          + PointT<T, sDim, tDim>(-z[1] * (1 - z[2]), (1 - z[0]) * (1 - z[2]), -(1 - z[0]) * z[1])
                * corners[3][k]
          + PointT<T, sDim, tDim>(-(1 - z[1]) * z[2], -(1 - z[0]) * z[2], (1 - z[0]) * (1 - z[1]))
                * corners[4][k]
          + PointT<T, sDim, tDim>((1 - z[1]) * z[2], -z[0] * z[2], z[0] * (1 - z[1]))
                * corners[5][k]
          + PointT<T, sDim, tDim>(z[1] * z[2], z[0] * z[2], z[0] * z[1]) * corners[6][k]
          + PointT<T, sDim, tDim>(-z[1] * z[2], (1 - z[0]) * z[2], (1 - z[0]) * z[1])
                * corners[7][k];
    }
    return TransformationT<T, sDim>(y);
  }
public:
  static const Point *LocalCornersHex();

  static const Point *LocalEdgesHex();

  static const Point *LocalFacesHex();

  static const Point *LocalFaceNormalsHex();

  static const Point *LocalCenterHex();

  static const vector<Rule> *Rhex;

  static const vector<Rule> *Rhex_barycentric;

  Hexahedron(const vector<Point> &z, int sd);

  virtual CELLTYPE Type() const override;

  CELLTYPE ReferenceType() const override;

  int Corners() const override;

  const Point &LocalCorner(int i) const override;

  const Point &LocalEdge(int i) const override;

  const Point &LocalFace(int i) const override;

  const Point &LocalFaceNormal(int i) const override;

  double LocalFaceArea(int i) const override;

  double FaceArea(int i) const override;

  const Point &LocalCenter() const override;

  virtual Point FaceLocalToGlobal(int face, const Point &) const override;

  virtual Point FaceLocalToLocal(int faceId, const Point &localOnFace) const override;

  virtual Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const override;

  virtual Point LocalToGlobal(const Point &z) const override;

  virtual Point GlobalToLocal(const Point &point) const override;

  Transformation GetTransformation(const Point &x) const override;

#ifdef BUILD_IA
  IAPoint GlobalToLocal(const IAPoint &point) const override;

  IAPoint LocalToGlobal(const IAPoint &z) const override;

  IAPoint FaceLocalToGlobal(int face, const IAPoint &) const override;

  IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const override;

  IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const override;

  IATransformation GetTransformation(const IAPoint &x) const override;

  IAInterval LocalFaceAreaIA(int faceId) const override;

  static const IAPoint *LocalFaceNormalsIAHex();

  const IAPoint &LocalFaceNormalIA(int faceId) const override;
#endif

  int Edges() const override;

  Point Edge(int i) const override;

  short edgecorner(unsigned short i, unsigned short j) const override;

  int Faces() const override;

  Point Face(int i) const override;

  CELLTYPE FaceType(int face) const override;

  short FarFace(int i) const override;

  short MiddleEdge(int i, int j) const override;

  short FaceCorners(int i) const override;

  short facecorner(unsigned short i, unsigned short j) const override;

  short FaceEdges(int i) const override;

  bool faceedgecircuit(int i) const override;

  short faceedge(int i, int j) const override;

  void GetRefinePoints(std::vector<Point> &z) const override;

  virtual const vector<Rule> &Refine(vector<Point> &z) const override;

  void GetRefineBarycentricPoints(std::vector<Point> &z, const Point &shiftCenter) const override;

  virtual const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                                     const Point &shiftCenter) const override;

  int dim() const override;

  bool plane() const override;

  static const Hexahedron *ReferenceHexahedron();

  const Cell *ReferenceCell() const override;

  int Children() const override;

  Point Child(int i) const override;

  bool PointInCell(const Point &global) const {
    Point local = GlobalToLocal(global);
    for (int i = 0; i < 3; i++) {
      if (local[i] < 0 || local[i] > 1) { return false; }
    }
    return true;
  }
};

// ----------------------------------------------------
//    HEXAHEDRON20
// ----------------------------------------------------

class Hexahedron20 : public Hexahedron {
public:
  static const vector<Rule> *R20hex;

  Hexahedron20(const vector<Point> &z, int sd);

  Point FaceLocalToGlobal(int face, const Point &) const override{
      THROW("Check if FaceLocalToGlobal needs to be implemented")}

  Point FaceLocalToLocal(int faceId, const Point &localOnFace) const override{
      THROW("Check if FaceLocalToLocal needs to be implemented")}

  Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const override{
      THROW("Check if FaceLocalToLocal needs to be implemented")}

  Point GlobalToLocal(const Point &z) const override{
      THROW("Check if GlobalToLocal needs to be implemented")}

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

  IATransformation GetTransformation(const IAPoint &x) const override{
      THROW("GetTransformation not implemented for IA types")}
#endif

  Point LocalToGlobal(const Point &z) const override;

  Transformation GetTransformation(const Point &x) const override;

  CELLTYPE Type() const override;

  void GetRefinePoints(std::vector<Point> &z) const override;

  const vector<Rule> &Refine(vector<Point> &z) const override;

  const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                             const Point &shiftCenter) const override {
    THROW("RefineBarycentric not implemented")
  }

  static const Hexahedron20 *ReferenceHexahedron20();
};

// ----------------------------------------------------
//    HEXAHEDRON27
// ----------------------------------------------------
class Hexahedron27 : public Hexahedron {
public:
  static const vector<Rule> *R27hex;

  Hexahedron27(const vector<Point> &z, int sd);

  Point FaceLocalToGlobal(int face, const Point &) const override{
      THROW("Check if FaceLocalToGlobal needs to be implemented")}

  Point FaceLocalToLocal(int faceId, const Point &localOnFace) const override{
      THROW("Check if FaceLocalToLocal needs to be implemented")}

  Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const override{
      THROW("Check if FaceLocalToLocal needs to be implemented")}

  Point GlobalToLocal(const Point &z) const override{
      THROW("Check if GlobalToLocal needs to be implemented")}

  Point LocalToGlobal(const Point &z) const override;

  Transformation GetTransformation(const Point &x) const override;

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

  IATransformation GetTransformation(const IAPoint &x) const override{
      THROW("GetTransformation not implemented for IA types")}
#endif

  CELLTYPE Type() const override;

  void GetRefinePoints(std::vector<Point> &z) const override;

  const vector<Rule> &Refine(vector<Point> &z) const override;

  const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                             const Point &shiftCenter) const override {
    THROW("RefineBarycentric not implemented")
  }

  static const Hexahedron27 *ReferenceHexahedron27();
};

#endif // HEXAHEDRON_H
