#ifndef TCELL_HPP
#define TCELL_HPP

#include <memory>
#include "Cell.hpp"
#include "Face.hpp"
#include "Intval.hpp"

class TimeInterval {
  double a;
  double b;
public:
  TimeInterval(double A, double B) : a(A), b(B) {}

  double Center() const { return 0.5 * (a + b); }

  double min() const { return a; }

  double max() const { return b; }
};

class TCell : public Cell {
  std::unique_ptr<const Cell> c;
  std::unique_ptr<const Intval> I;
private:
  vector<Point> createCorners(const Cell *cell, const TimeInterval interval) {
    vector<Point> corners(2 * cell->Corners());
    for (int i = 0; i < cell->Corners(); i++) {
      corners[i] = cell->Corner(i).CopyWithT(interval.min());
      corners[i + cell->Corners()] = cell->Corner(i).CopyWithT(interval.max());
    }
    return corners;
  }
public:
  TCell(const Cell *cell, const TimeInterval &interval) :
      Cell(createCorners(cell, interval), 0), c(std::unique_ptr<const Cell>(cell)) {
    std::vector<Point> pts{interval.min(), interval.max()};
    I = std::make_unique<const Intval>(pts, 0);
  }

  void GetRefinePoints(vector<Point> &z) const override { c->GetRefinePoints(z); }

  const vector<Rule> &Refine(vector<Point> &z) const override { return c->Refine(z); }

  int dim() const override { return c->dim(); }

  bool plane() const override { return true; }

  CELLTYPE Type() const override { return SpaceTimeCellType(c->Type()); }

  CELLTYPE ReferenceType() const override { return SpaceTimeCellType(c->ReferenceType()); }

  const Cell *ReferenceCell() const override { return c->ReferenceCell(); }

  const Point &LocalCenter() const override { THROW("TCell: LocalCenter not implemented.") }

  int Corners() const override { return corners.size(); }

  const Point &LocalCorner(int i) const override {
    THROW("TCell::LocalCorner not implemented.");
    return c->LocalCorner(i);
  }

  int Edges() const override {
    return 2 * c->Edges() + c->Corners();
    Exit("not implemented.");
  }

  Point Edge(int i) const override {
    if (i < c->Edges()) {
      return c->Edge(i).CopyWithT(min());
    } else if (i < c->Edges() + c->Corners()) {
      return c->Corner(i - c->Edges()).CopyWithT(I->Center()[0]);
    } else {
      return c->Edge(i - c->Corners() - c->Edges()).CopyWithT(max());
    }
  }

  const Point &LocalEdge(int) const override { THROW("TCell::LocalEdge not implemented."); }

  short edgecorner(unsigned short i, unsigned short ci) const override {
    if (i < c->Edges()) {
      return c->edgecorner(i, ci);
    } else if (i < c->Edges() + c->Corners()) {
      return i - c->Edges() + ci * c->Corners();
    } else {
      return c->edgecorner(i - c->Corners() - c->Edges(), ci) + c->Corners();
    }
  }

  int Faces() const override { return c->Faces() + 2; }

  Point Face(int i) const override {
    if (i < Faces() - 2) return Point(c->Face(i), I->Center()[0]);
    else if (i == Faces() - 1) return Point(Center(), max());
    else return Point(Center(), min());
  }

  const Point &LocalFace(int i) const override {
    THROW("TCell::LocalFace for TCell not implemented")
    return c->LocalFace(i);
  }

  CELLTYPE FaceType(int face) const override { THROW("TCell::FaceType for TCell not implemented") }

  const Point &LocalFaceNormal(int i) const override {
    THROW("TCell::LocalFaceNormal not implemented");
    return c->LocalFaceNormal(i);
  }

  double LocalFaceArea(int face) const override {
    THROW("TCell::LocalFaceArea not implemented");
    return c->LocalFaceArea(face);
  }

  double FaceArea(int face) const override {
    THROW("TCell::FaceArea not implemented");
    return c->FaceArea(face);
  }

  short FaceCorners(int f) const override {
    if (c->dim() == 3) { return 8; }
    if (f < c->Faces()) {
      // Todo 3D
      return 4;
    } else {
      return c->Corners();
    }
  }

  short facecorner(unsigned short f, unsigned short c_index) const override {
    if (f < c->Faces()) {
      switch (c_index) {
      case 0:
        return f;
      case 1:
        return f + c->Corners();
      case 2:
        return (f + 1) % c->Faces() + c->Corners();
      case 3:
        return (f + 1) % c->Faces();
      default:
        THROW("no");
      }
    } else if (f == c->Faces()) {
      return c_index;
    } else {
      return c_index + c->Corners();
    }
  }

  short FaceEdges(int f) const override {
    return 0;
    if (f < c->Faces()) {
      // Todo 3D
      return 4;
    } else {
      return c->Faces();
    }
  }

  short faceedge(int f, int e_index) const override {
    if (f < c->Faces()) {
      switch (e_index) {
      case 0:
        return c->Corners() + f;
      case 1:
        return f;
      case 2:
        return c->Corners() + c->Faces() + f;
      case 3:
        return (f + 1) % c->Corners();
      default:
        THROW("no");
      }
    } else if (f == c->Faces()) {
      return e_index + c->Faces();
    } else {
      return e_index + 2 * c->Faces();
    }
  }

  int Children() const override {
    THROW("TCell::Children not implemented.")
    return c->Children();
  }

  Point Child(int i) const override {
    THROW("TCell::Child not implemented.")
    return c->Child(i);
  }

  Point LocalToGlobal(const Point &localPoint) const override {
    Point spaceGlobal = c->LocalToGlobal(localPoint);
    Point timeGlobal = I->LocalToGlobal(Point(localPoint.t()));
    return spaceGlobal.WithT(timeGlobal[0]);
  }

  Point GlobalToLocal(const Point &globalPoint) const override {
    Point spaceLocal = c->GlobalToLocal(globalPoint);
    Point timeLocal = I->GlobalToLocal(Point(globalPoint.t()));
    return spaceLocal.WithT(timeLocal[0]);
  }

  Point FaceLocalToGlobal(int face, const Point &) const override{
      THROW("TCell::FaceLocalToGlobal not implemented.")}

  Point FaceLocalToLocal(int face, const Point &) const override{
      THROW("TCell::FaceLocalToLocal not implemented.")}

  Point FaceLocalToLocalOriented(int face, const Point &) const override {
    THROW("TCell::FaceLocalToLocalOriented not implemented.")
  }

  virtual bool IsSpaceTime() const { return true; }

  double min() const override { return I->Corner(0)[0]; }

  double max() const override { return I->Corner(1)[0]; }

  const Cell &SpaceCell() const override { return *c; }

  const Cell &TimeCell() const override { return *I; }

  std::vector<Point> getChildrenInSpace() const {
    vector<Point> Children_in_space(Children());
    for (int l = 0; l < Children(); ++l) {
      Children_in_space[l] = Child(l).CopyWithT((*this)().t());
    }
    return Children_in_space;
  }

  std::vector<Point> getChildrenInTime() const {
    return {Center().CopyWithT(0.5 * (min() + Center().t())),
            Center().CopyWithT(0.5 * (Center().t() + max()))};
  }

  Transformation GetTransformation(const Point &z) const override {
    auto SpaceT = c->GetTransformation(z);
    double time_len = max() - min();
    return {SpaceT, time_len};
  }

  virtual bool PointInCell(const Point &global) const {
    bool in_space = c->PointInCell(global.CopyWithT(0));
    bool in_time = I->PointInCell({global.t()});
    return in_space && in_time;
  }

  ~TCell() override = default;
};

inline std::ostream &operator<<(std::ostream &os, const TCell &C) { return C.print(os); }

inline std::ostream &operator<<(std::ostream &os, const TCell *C) { return os << *C; }

#endif // TCELL_HPP
