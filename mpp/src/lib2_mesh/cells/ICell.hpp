#ifndef ICELL_H
#define ICELL_H

#include "Rule.hpp"
#include "Transformation.hpp"

#ifdef USE_DATAMESH

#include "DataSet.hpp"

#endif

template<typename T, int sDim, int tDim>
PointT<T, sDim, tDim> calculateCenter(const std::vector<PointT<T, sDim, tDim>> &corners) {
  PointT<T, sDim, tDim> center = Origin;
  for (const auto &corner : corners) {
    center += corner;
  }
  center /= corners.size();
  return center;
}

/**
 * \class Cell
 * \brief Stores the geometry of a single cell
 *
 * \param corners     The corners of the cell as read from the geometry file
 * \param subdomain   The subdomain of the cell as read from the geometry file
 * \param midPoint    The midpoint of the cell
 *
 * Cells defined in the geometry file are stored as Cell objects in the Mesh. Each Cell stores the
 * location of its corners, edge-centers, face-centers and center point. When inheriting from Cell,
 * is is crucial to make sure the orientation of the created cell objects are always the same. As
 * the corners obtained from the geometry don't have to be oriented correctly, a Cell should first
 * reorder its corners (if necessary) and then calculated faces, ...
 */
class Cell {
  /// Stores problem-specific data, such as encoded material parameters, densities, ...
  short subdomain;

  friend class Mesh;
protected:
  std::vector<Point> corners;
  Point midPoint;

#ifdef USE_DATAMESH
  DataContainer data;
#endif

  Cell(std::vector<Point> x, short sd) :
      subdomain(sd), corners(std::move(x)), midPoint(calculateCenter(corners)) {}
public:
  virtual ~Cell() {}

  /// Fills the vector with the points needed for refinement of this cell
  virtual void GetRefinePoints(std::vector<Point> &) const = 0;

  /// Returns a vector containing the rules for the refinement procedure of this cell
  virtual const std::vector<Rule> &Refine(std::vector<Point> &) const = 0;

  /// Fills the vector with the points needed for refinement of this cell
  virtual void GetRefineBarycentricPoints(std::vector<Point> &, const Point &shiftCenter) const {
    THROW("Barycentric refinement not implemented for cell type " + std::to_string(Type()))
  }

  /// Returns a vector containing the rules for barycentric refinement procedure of this cell
  virtual const std::vector<Rule> &RefineBarycentric(std::vector<Point> &,
                                                     const Point &shiftCenter) const {
      THROW("Barycentric refinement not implemented for cell type " + std::to_string(Type()))}

  /// Returns a copy of the corners if this cell
  std::vector<Point> copyCorners() const {
    return std::vector<Point>(corners);
  }

  /// Returns a reference to the corners of this cell
  const std::vector<Point> &AsVector() const { return corners; }

  /// Returns the number of corners contained in corners. Note that this might differs from the
  /// number returned by the function Corners() (e.g. if one uses tetrahedrons including edges
  /// (T10 elements))
  int Size() const { return corners.size(); }

  /// Returns the dimension of this cell
  virtual int dim() const = 0;

  /// Returns true, if the cell is planar
  virtual bool plane() const = 0;

  /// Returns the cell type
  virtual CELLTYPE Type() const = 0;

  /// Returns the basic cell type, i.e., if one uses tetrahedrons including edges (T10 elements),
  /// the basic cell type is still TETRAHEDRON
  virtual CELLTYPE ReferenceType() const = 0;

  /// Returns the reference cell for this cell type.
  virtual const Cell *ReferenceCell() const = 0;

  /// Returns the subdomain
  virtual short Subdomain() const { return subdomain; }

  /// Returns the Volume of the cell
  virtual double GetVolume() const {THROW("Cell::GetVolume not implemented")};

  /// Returns the global center point of the cell
  Point Center() const { return midPoint; }

  /// Returns the global center point of the cell
  Point operator()() const { return midPoint; }

  /// Returns the local center point of the cell, i.e., of the reference cell
  virtual const Point &LocalCenter() const = 0;

  /// Returns the number of corners of this cell
  virtual int Corners() const = 0;

  /// Returns the cornerId-th global corner
  virtual const Point &Corner(int cornerId) const {
#if DebugLevel > 0
    try {
      return corners.at(cornerId);
    } catch (const std::out_of_range &oor) {
      THROW("Corner " + to_string(cornerId) + " does not exist.")
    }
#else
    return corners[cornerId];
#endif
  }

  /// Returns the cornerId-th global corner
  virtual const Point &operator[](int cornerId) const { return Corner(cornerId); }

  /// Returns the cornerId-th local corner, i.e., of the reference cell
  virtual const Point &LocalCorner(int cornerId) const = 0;

  /// Returns the number of edges of this cell
  virtual int Edges() const = 0;

  /// Returns the global midpoint of the edgeId-th edge
  virtual Point Edge(int edgeId) const = 0;

  /// Returns the local midpoint of the edgeId-th edge, i.e., of the reference cell
  virtual const Point &LocalEdge(int edgeId) const = 0;

  /// Returns the amount of corners of the edgeId-th edge (at most 2, only different for Intval)
  virtual int EdgeCorners(int edgeId) const { return 2; }

  /// Returns the cornerId of the edgecornerId-th edge corner of the edgeId-th edge
  virtual short edgecorner(unsigned short edgeId, unsigned short edgecornerId) const = 0;

  /// Returns the global edgecornerId-th edge corner of the edgeId-th edge
  const Point &EdgeCorner(unsigned short edgeId, unsigned short edgecornerId) const {
    return Corner(edgecorner(edgeId, edgecornerId));
  };

  /// Returns the number of faces of this cell
  virtual int Faces() const = 0;

  /// Returns the global midpoint of the faceId-th face
  virtual Point Face(int faceId) const = 0;

  /// Returns the local midpoint of the faceId-th face, i.e., of the reference cell
  virtual const Point &LocalFace(int faceId) const = 0;

  /// Return the faceId of the face corresponding to the faceMidpoint
  int Face(const Point &faceMidpoint) const {
    for (int faceId = 0; faceId < Faces(); ++faceId)
      if (Face(faceId) == faceMidpoint) return faceId;
    THROW("Face not found")
  }

  /// Returns the cell type of the faceId-th face
  virtual CELLTYPE FaceType(int faceId) const = 0;

  /// Return the local face normal of the faceId-th face, i.e., on the reference cell
  virtual const Point &LocalFaceNormal(int faceId) const = 0;

  template<typename T, int sDim, int tDim>
  const PointT<T, sDim, tDim> &LocalFaceNormalT(int faceId) const;

  /// Return the face area of the global faceId-th face
  virtual double FaceArea(int faceId) const = 0;

  /// Return the face area of the local faceId-th face, i.e., of the reference cell
  virtual double LocalFaceArea(int faceId) const = 0;

  template<typename T>
  T LocalFaceAreaT(int faceId) const;

  virtual VectorField GetNablaNi(int faceId) const {THROW("Cell::GetNablaNi not implemented")};

  /// Returns the amount of corners of the faceId-th face
  virtual short FaceCorners(int faceId) const = 0;

  /// Returns the cornerId of the facecornerId-th face corner of the faceId-th face
  virtual short facecorner(unsigned short faceId, unsigned short facecornerId) const = 0;

  /// Returns the global facecornerId-th face corner of the faceId-th face
  const Point &FaceCorner(unsigned short faceId, unsigned short facecornerId) const {
    return Corner(facecorner(faceId, facecornerId));
  };

  /// Returns the local facecornerId-th face corner of the faceId-th face
  const Point &LocalFaceCorner(unsigned short faceId, unsigned short facecornerId) const {
    return LocalCorner(facecorner(faceId, facecornerId));
  };

  /// Returns the amount of edges of the faceId-th face
  virtual short FaceEdges(int faceId) const = 0;

  /// Returns the edgeId of the faceedgeId-th face edge of the faceId-th face
  virtual short faceedge(int faceId, int faceedgeId) const = 0;

  /// Returns the global midpoint of the faceedgeId-th edge of the faceId-th face
  virtual Point FaceEdge(int faceId, int faceedgeId) const {
    return Edge(faceedge(faceId, faceedgeId));
  }

  /// Returns the local midpoint of the faceedgeId-th edge of the faceId-th face
  virtual Point LocalFaceEdge(int faceId, int faceedgeId) const {
    return LocalEdge(faceedge(faceId, faceedgeId));
  }

  /// Returns the cornerId of the faceedgecornerId-th edge of the faceedgeId-th face edge
  /// of the faceId-th face
  short faceedgecorner(int faceId, int faceedgeId, int faceedgecornerId) const {
    return edgecorner(faceedge(faceId, faceedgeId), faceedgecornerId);
  }

  /// Returns the global faceedgecornerId-th edge corner of the faceedgeId-th face edge
  /// of the faceId-th face
  const Point &FaceEdgeCorner(int faceId, int faceedgeId, int faceedgecornerId) const {
    return EdgeCorner(faceedge(faceId, faceedgeId), faceedgecornerId);
  }

  /// Returns the amount of children of this cell
  virtual int Children() const = 0;

  /// Returns the global center point of the childId-th child
  virtual Point Child(int childId) const = 0;

  /// Returns the transformation at the given point global
  virtual Transformation GetTransformation(const Point &global = Origin) const = 0;

  /// Transforms a given local point of the reference cell to a global point in the cell
  virtual Point LocalToGlobal(const Point &local) const = 0;

  /// Transforms a given local point of the reference cell to a global point in the cell
  Point operator[](const Point &local) const { return LocalToGlobal(local); }

  /// Transforms a given global point of the cell to a local point in the reference cell
  virtual Point GlobalToLocal(const Point &global) const = 0;

  /// Transforms a given local point of the faceId-th face (in local face coordinates)
  /// to a global point in the cell
  virtual Point FaceLocalToGlobal(int faceId, const Point &localOnFace) const = 0;

  /// Transforms a given local point of the faceId-th face (in local face coordinates)
  /// to a local point in the cell
  virtual Point FaceLocalToLocal(int faceId, const Point &localOnFace) const = 0;

  /// Transforms a given local point of the faceId-th face (in local face coordinates)
  /// to a local point in the cell such that on neighbouring cells the points coincide
  virtual Point FaceLocalToLocalOriented(int faceId, const Point &localOnFace) const = 0;

  /// Checks if a given global point is contained in the cell
  virtual bool PointInCell(const Point &global) const { THROW("Cell::PointInCell not implemented") }

  /**
   * Interval arithmetic functions (can be removed if cells are template based)
   */
#ifdef BUILD_IA
  /// Returns the ia transformation at the given point global
  virtual IATransformation GetTransformation(const IAPoint &global) const {
    THROW("GetTransformation not implemented for interval arithmetic")
  }

  /// Transforms a given global ia point of the cell to a local point in the reference cell
  virtual IAPoint GlobalToLocal(const IAPoint &global) const {
    THROW("GlobalToLocal not implemented for interval arithmetic")
  }

  /// Transforms a given local ia point of the reference cell to a global ia point in the cell
  virtual IAPoint LocalToGlobal(const IAPoint &local) const {
    THROW("LocalToGlobal not implemented for interval arithmetic")
  }

  /// Transforms a given local ia point of the faceId-th face (in local face coordinates)
  /// to a global ia point in the cell
  virtual IAPoint FaceLocalToGlobal(int faceId, const IAPoint &localOnFace) const {
    THROW("FaceLocalToGlobal not implemented for interval arithmetic")
  }

  /// Transforms a given local ia point of the faceId-th face (in local face coordinates)
  /// to a local ia point in the cell
  virtual IAPoint FaceLocalToLocal(int faceId, const IAPoint &localOnFace) const {
    THROW("FaceLocalToLocal not implemented for interval arithmetic")
  }

  /// Transforms a given local ia point of the faceId-th face (in local face coordinates)
  /// to a local ia point in the cell such that on neighbouring cells the points coincide
  virtual IAPoint FaceLocalToLocalOriented(int faceId, const IAPoint &localOnFace) const {
    THROW("FaceLocalToLocalOriented not implemented for interval arithmetic")
  }

  /// Return the ia face area of the local faceId-th face, i.e., of the reference cell
  virtual IAInterval LocalFaceAreaIA(int faceId) const {
    THROW("LocalFaceAreaIA not implemented for interval arithmetic")
  }

  /// Return the ia local face normal of the faceId-th face, i.e., on the reference cell
  virtual const IAPoint &LocalFaceNormalIA(int faceId) const {
    THROW("LocalFaceNormalIA not implemented for interval arithmetic")
  }
#endif

  /**
   * Space time functions
   */
  virtual bool IsSpaceTime() const { return false; }

  virtual double min() const { return 0.0; }

  virtual double max() const { return 0.0; }

  virtual const Cell &SpaceCell() const { return *this; }

  virtual const Cell &TimeCell() const { THROW("Cell::TimeCell not impl.") }

#ifdef USE_SPACETIME
  virtual std::vector<Point> getChildrenInSpace() const { THROW("Cell::ChildrenInSpace not impl.") }

  virtual std::vector<Point> getChildrenInTime() const { THROW("Cell::ChildrenInSpace not impl.") }

  /*virtual Transformation GetSpaceTimeTransformation(const Point &q = Origin) const {
    THROW("Cell::GetSpaceTimeTransformation not implemented.")
  };*/
#endif

  /**
   * Data mesh functions
   */
#ifdef USE_DATAMESH
  const DataContainer &GetData() const { return data; }

  void SetData(const DataContainer &d) { data = d; }
#endif

  template<typename OSTREAM>
  OSTREAM &print(OSTREAM &os) const {
    for (int i = 0; i < Size(); ++i)
      os << (*this)[i] << "|";
    return os << " tp " << int((*this).Type()) << " sd " << (*this).Subdomain();
  }

  /**
   * Not sure what the following functions actually do
   */
  virtual short FarFace(int) const { return -1; }

  virtual short MiddleEdge(int i, int j) const { return -1; }

  virtual bool faceedgecircuit(int i) const { return false; }
};

inline std::ostream &operator<<(std::ostream &os, const Cell &C) { return C.print(os); }

inline std::ostream &operator<<(std::ostream &os, const Cell *C) { return os << *C; }

inline const Point *ReferenceEdges(const Cell &C) {
  Point *z = new Point[C.Edges()];
  for (int i = 0; i < C.Edges(); ++i)
    z[i] = C.LocalToGlobal(
        0.5 * (C.LocalCorner(C.edgecorner(i, 0)) + C.LocalCorner(C.edgecorner(i, 1))));
  return z;
}

inline const Point *ReferenceFaces(const Cell &C) {
  Point *z = new Point[C.Faces()];
  for (int i = 0; i < C.Faces(); ++i) {
    z[i] = Origin;
    for (int j = 0; j < C.FaceCorners(i); ++j) {
      z[i] += C.FaceCorner(i, j);
    }
    z[i] *= 1.0 / C.FaceCorners(i);
  }
  return z;
}

template<typename POINT>
inline const POINT *ReferenceFaceNormals(const Cell &C) {
  POINT *z = new POINT[C.Faces()];
  for (int i = 0; i < C.Faces(); ++i) {
    z[i] = Origin;
    if (C.FaceCorners(i) == 2) {
      const Point &P0 = C.FaceCorner(i, 0);
      const Point &P1 = C.FaceCorner(i, 1);
      POINT N(P1[1] - P0[1], P0[0] - P1[0]);
      z[i] = (1.0 / norm(N)) * N;
    } else if (C.FaceCorners(i) > 2) {
      const Point &P0 = C.FaceCorner(i, 0);
      const Point &P1 = C.FaceCorner(i, 1);
      const Point &P2 = C.FaceCorner(i, 2);
      POINT X = P1 - P0;
      POINT Y = P2 - P0;
      POINT N = curl(X, Y);
      z[i] = (1.0 / norm(N)) * N;
    }
  }
  return z;
}

inline const Point *ReferenceCenter(const Cell &C) {
  Point *z = new Point;
  *z = Origin;
  for (int i = 0; i < C.Corners(); ++i)
    *z += C[i];
  *z *= (1.0 / C.Corners());
  return z;
}

#endif // ICELL_H
