#ifndef _CELL_H_
#define _CELL_H_

#include "ICell.hpp"

#include "DataSet.hpp"


Cell *CreateCell(CELLTYPE tp, int subdomain, const vector<Point> &corners, double a = 0.0,
                 double b = 0.0);

const Cell *ReferenceCell(CELLTYPE tp);

class cell : public std::unordered_map<Point, Cell *>::const_iterator {
  typedef std::unordered_map<Point, Cell *>::const_iterator Iterator;
public:
  cell() {}

  cell(Iterator c) : Iterator(c) {}

  const vector<Point> &AsVector() { return (*this)->second->AsVector(); }

  int size() const { return (*this)->second->Size(); }

  const Point &operator()() const { return (*this)->first; }

  const Point &Center() const { return (*this)->first; }

  const Cell &operator*() const { return *((*this)->second); }

  CELLTYPE Type() const { return (*this)->second->Type(); }

  CELLTYPE ReferenceType() const { return (*this)->second->ReferenceType(); }

  const Cell *ReferenceCell() { return (*this)->second->ReferenceCell(); }

  short Subdomain() const { return (*this)->second->Subdomain(); }

  double GetVolume() const { return (*this)->second->GetVolume(); }

  VectorField GetNablaNi(int i) const { return (*this)->second->GetNablaNi(i); }

  int Corners() const { return (*this)->second->Corners(); }

  const Point &Corner(int i) const { return (*this)->second->Corner(i); }

  const Point &LocalCorner(int i) const { return (*this)->second->LocalCorner(i); }

  const Point &LocalEdge(int i) const { return (*this)->second->LocalEdge(i); }

  const Point &LocalFace(int i) const { return (*this)->second->LocalFace(i); }

  const Point &LocalCenter() const { return (*this)->second->LocalCenter(); }

  const Point &LocalFaceNormal(int i) const { return (*this)->second->LocalFaceNormal(i); }

  Point LocalToGlobal(const Point &z) const { return (*this)->second->LocalToGlobal(z); }

  Point GlobalToLocal(const Point &z) const { return (*this)->second->GlobalToLocal(z); }

  bool PointInCell(const Point &z) const { return (*this)->second->PointInCell(z); }

  Point operator[](const Point &z) const { return (*this)->second->LocalToGlobal(z); }

  Transformation GetTransformation(const Point &z = Origin) const {
    return (*this)->second->GetTransformation(z);
  }


#ifdef BUILD_IA
  IAPoint LocalToGlobal(const IAPoint &z) const { return (*this)->second->LocalToGlobal(z); }

  IATransformation GetTransformation(const IAPoint &z) const {
    return (*this)->second->GetTransformation(z);
  }
#endif

  virtual bool IsSpaceTime() const { return (*this)->second->IsSpaceTime(); }

  double min() const { return (*this)->second->min(); }

  double max() const { return (*this)->second->max(); }

  std::vector<Point> copyCorners() const { return (*this)->second->copyCorners(); }

  virtual const Cell &SpaceCell() const { return (*this)->second->SpaceCell(); }

  virtual const Cell &TimeCell() const { return (*this)->second->TimeCell(); }

#ifdef USE_SPACETIME
  std::vector<Point> getChildrenInSpace() const { return (*this)->second->getChildrenInSpace(); }

  std::vector<Point> getChildrenInTime() const { return (*this)->second->getChildrenInTime(); }
#endif

  const Point &operator[](int i) const { return (*this)->second->Corner(i); }

  int Edges() const { return (*this)->second->Edges(); }

  Point Edge(int i) const { return (*this)->second->Edge(i); }

  const Point &EdgeCorner(int i, int j) const { return (*this)->second->EdgeCorner(i, j); }

  short edgecorner(int i, int j) const { return (*this)->second->edgecorner(i, j); }

  int Faces() const { return (*this)->second->Faces(); }

  Point Face(int i) const { return (*this)->second->Face(i); }

  short FarFace(int i) const { return (*this)->second->FarFace(i); }

  short MiddleEdge(int i, int j) const { return (*this)->second->MiddleEdge(i, j); }

  bool faceedgecircuit(int i) const { return (*this)->second->faceedgecircuit(i); }

  short FaceCorners(int i) const { return (*this)->second->FaceCorners(i); }

  const Point &FaceCorner(int i, int j) const { return (*this)->second->FaceCorner(i, j); }

  short facecorner(int i, int j) const { return (*this)->second->facecorner(i, j); }

  short FaceEdges(int i) const { return (*this)->second->FaceEdges(i); }

  const Point &FaceEdgeCorner(int i, int j, int k) const {
    return (*this)->second->FaceEdgeCorner(i, j, k);
  }

  Point FaceEdge(int i, int j) const { return (*this)->second->FaceEdge(i, j); }

  short faceedge(int i, int j) const { return (*this)->second->faceedge(i, j); }

  short faceedgecorner(int i, int j, int k) const {
    return (*this)->second->faceedgecorner(i, j, k);
  }

  short facecorner(const Point &z) const {
    for (int i = 0; i < Faces(); ++i)
      if (z == Face(i)) return i;
    return -1;
  }

  int FaceId(const Point &z) const {
    for (int i = 0; i < Faces(); ++i)
      if (Face(i) == z) return i;
    THROW("Faceid should not get here.")
  }

  double LocalFaceArea(int face) const { return (*this)->second->LocalFaceArea(face); }

  double FaceArea(int face) const { return (*this)->second->FaceArea(face); }

  const vector<Rule> &Refine(vector<Point> &z) const { return (*this)->second->Refine(z); }

  const vector<Rule> &RefineBarycentric(vector<Point> &z, const Point &shiftCenter) const {
    return (*this)->second->RefineBarycentric(z, shiftCenter);
  }

  int dim() const { return (*this)->second->dim(); }

  bool plane() const { return (*this)->second->plane(); }

  const Point &LocalEdgeCorner(int i, int k) const {
    return LocalCorner((*this)->second->edgecorner(i, k));
  }

  const Point &LocalFaceCorner(int i, int k) const {
    return LocalCorner((*this)->second->facecorner(i, k));
  }

  int Children() const { return (*this)->second->Children(); }

  Point Child(int i) const { return (*this)->second->Child(i); }

  template<typename OSTREAM>
  OSTREAM &print(OSTREAM &os) const {
    return os << (*this)->first << " : " << *((*this)->second);
  }

#ifdef USE_DATAMESH

  const DataContainer &GetData() const { return (*this)->second->GetData(); }

  void SetData(const DataContainer &container) { (*this)->second->SetData(container); }

#endif

  bool operator<(const cell &other) const { return (*this)->first < other->first; }
};

class Buffer;

Buffer &operator<<(Buffer &b, const Cell &c);

Buffer &operator<<(Buffer &b, const cell &c);

struct GetCell {};

Buffer &operator>>(Buffer &b, Cell *&c);

inline std::ostream &operator<<(std::ostream &os, const cell &c) { return c.print(os); }

inline std::ostream &operator<<(std::ostream &s, const std::vector<std::list<cell>> &t) {
  for (int n = 0; n < t.size(); n++) {
    s << n << ":\n";
    for (const cell &c : t[n]) {
      s << c() << "\n";
    }
  }
  s << std::endl;
  return s;
}

bool LessCell_x(const cell &c0, const cell &c1);

bool LessCell_y(const cell &c0, const cell &c1);

bool LessCell_z(const cell &c0, const cell &c1);

#endif
