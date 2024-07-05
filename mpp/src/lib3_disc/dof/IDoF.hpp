#ifndef _IDOF_H_
#define _IDOF_H_

#include <map>
#include <string>
#include "Mesh.hpp"

struct DegreePair {
  short space;
  short time;

  DegreePair decreaseBoth() const {
    return {(short)(MAX(space - 1, 0)), (short)(MAX(time - 1, 0))};
  }

  DegreePair increaseBoth() const {
    return {(short)(MAX(space - 1, 0)), (short)(MAX(time - 1, 0))};
  }
};

inline bool operator==(const DegreePair &lhs, const DegreePair &rhs) {
  return lhs.space == rhs.space && lhs.time == rhs.time;
}

using DegreePairProduct = std::pair<DegreePair, DegreePair>;

namespace std {
template<>
struct hash<DegreePair> {
  inline size_t operator()(const DegreePair &x) const { return x.space * 10 + x.time; }
};

template<>
struct hash<DegreePairProduct> {
  inline size_t operator()(const DegreePairProduct &x) const {
    return hash<DegreePair>()(x.first) * 10 + hash<DegreePair>()(x.second);
  }
};

} // namespace std

// TODO: What are different PointTypes???
constexpr int MaxPointTypes = 10;

enum NODETYPE { VERTEX = 0, EDGE = 1, FACE = 2, CELL = 3 };

class SubVectorMask {
  int tp[MaxPointTypes];
  bool *m[MaxPointTypes];
public:
  SubVectorMask();

  SubVectorMask(const int *n, const char *c);

  SubVectorMask(const SubVectorMask &s);

  SubVectorMask &operator=(const SubVectorMask &s);

  ~SubVectorMask();

  const bool *operator[](int i) const { return m[i]; }

  friend std::ostream &operator<<(std::ostream &s, const SubVectorMask &S);
};

class IDoF {
  int _n_infty;
  std::map<std::string, SubVectorMask> sub;
  bool bnd;
public:
  IDoF(int n = 0, bool b = false) : _n_infty(n), bnd(b) {}

  virtual ~IDoF() {}

  virtual std::string Name() const { THROW("Name() not implemented") }

  virtual std::string Name(const Cell &) const { THROW("Name(const Cell&) not implemented") }

  /// Override this function e.g. in MixedDoF
  virtual short NumberOfDoFs() const { return 1; }

  virtual short NumberOfComponents() const { THROW("NumberOfComponents() not implemented") }

  /// Returns the true nodal points and NOT the storage points
  virtual short NumberOfNodalPoints(const Cell &c) const = 0;

  virtual short NumberOfNodalPoints(const int sDeg, const int tDeg) const {
    THROW("NumberOfNodalPoints")
  }

  virtual std::vector<Point> GetNodalPoints(const Cell &c) const = 0;

  virtual std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const = 0;

  virtual short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const = 0;

  virtual short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const = 0;

  std::vector<Point> GetNodalPointsOnFace(const Cell &c, int faceId) const;

  virtual short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const = 0;

  virtual short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const = 0;

  std::vector<Point> GetNodalPointsOnEdge(const Cell &c, int edgeId) const;

  /// Returns the number of storage points required for setting up the MatrixGraph
  virtual short NumberOfStoragePoints(const Cell &c) const;

  /// Returns the storage points required for setting up the MatrixGraph
  virtual std::vector<Point> GetStoragePoints(const Cell &c) const;

  /// Returns the amount of memory needed at the storage points
  virtual std::vector<short> AllocationSizesAtStoragePoints(const Cell &c) const;

  virtual std::vector<short> AccumulatedAllocationSizes(const Cell &c) const;

  virtual short NumberOfStoragePointsOnFace(const Cell &c, int faceId) const;

  virtual short IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const;

  std::vector<Point> GetStoragePointsOnFace(const Cell &c, int faceId) const;

  virtual short NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const;

  virtual short IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const;

  std::vector<Point> GetStoragePointsOnEdge(const Cell &c, int edgeId) const;

  // TODO: Maybe removable
  virtual void NodalPointsLocal(const int &, vector<Point> &, const int k = 0) const {
    THROW("Not Implemented")
  }

  // TODO: Maybe removable
  virtual int NodalPointsLocal(const int &, const int k = 0) const { THROW("Not Implemented") }

  virtual void SetDegree(Point p, DegreePair pair) { THROW("SetDegree") }

  virtual DegreePair GetDegree(Point p) const { THROW("GetDegree") }

  virtual DegreePair GetDegree(const Cell &) const { THROW("GetDegree") }

  virtual DegreePair GetMaxDegree() const { THROW("GetMaxDegree") }

  virtual int get_space_deg(const Point &) const { THROW("Not Implemented") }

  virtual int get_cell_deg(const Cell &) const { THROW("Not Implemented") }

  virtual int get_time_deg(const Point &) const { THROW("Not Implemented") }

  virtual void set_space_deg(const Point &, int) {}

  virtual void set_cell_deg(const Cell &, int) {}

  virtual void set_time_deg(const Point &, int) {}

  virtual int get_dim() const { return -1; }

  virtual void increase_cell_deg(const Cell &) {}

  virtual void decrease_cell_deg(const Cell &) {}

  virtual int get_time_deg(const Cell &) const { THROW("Not Implemented") }

  virtual int get_space_deg(const Cell &) const { THROW("Not Implemented") }

  virtual void set_space_deg(const Cell &, int) {}

  virtual void set_time_deg(const Cell &, int) {}

  virtual void increase_space_deg(const Cell &) {}

  virtual void decrease_space_deg(const Cell &) {}

  virtual void increase_time_deg(const Cell &) {}

  virtual void decrease_time_deg(const Cell &) {}

  virtual void communicate() {}

  virtual int get_m(const Point &) const { return 0; }

  virtual int TypeDoF(int i) const { THROW("Not implemented") }

  void TypeDoFs(int *n) const;

  int nsub() const { return sub.size(); }

  int ninf() const { return _n_infty; }

  int n_infty() const { return _n_infty; }

  bool Bnd() const { return bnd; }

  const std::map<std::string, SubVectorMask> &Sub() const { return sub; }

  const SubVectorMask &GetSubVector(const char *name) const;

  void AddSubVector(const char *name, const char *c, const int *n);

  void AddSubVector(const char *name, const char *c);

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const IDoF &D);
};

#endif // of #ifndef _IDOF_H_
