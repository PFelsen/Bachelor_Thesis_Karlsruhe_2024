#ifndef SPACETIMEMATRIXGRAPH_HPP
#define SPACETIMEMATRIXGRAPH_HPP

#include "IMatrixGraph.hpp"

class SpaceTimeMatrixGraph : public IMatrixGraph {
protected:
  bool isBlockDiagonal = false;

  bool isDgInTime;

  DegreePair defaultDegree;

  DegreePair max_deg;

  int LOWEST_POSSIBLE_TIME_DEG = 0;

  std::unordered_map<Point, DegreePair> cell_deg;

  void AddNBCells();

  void InitProcsets();

  void Init() override;

  void STOverlapCellsWithCorners();

  void STOverlapCellsWithFaces();
public:
  void communicate();

  void SetDegree(Point p, DegreePair pair) {
    cell_deg[p] = pair;
    if (pair.space > max_deg.space) max_deg.space = pair.space;
    if (pair.time > max_deg.time) max_deg.time = pair.time;
  }

  DegreePair GetDegree(Point p) const { return cell_deg.find(p)->second; }

  DegreePair GetDegree(const Cell &c) const { return GetDegree(c()); }

  DegreePair GetMaxDegree() const { return max_deg; }

  string Name(const Cell &c) const {
    DegreePair degs = cell_deg.at(c());
    std::stringstream ss;
    ss << "AdaptiveSTDoF Base sd=" << degs.space << " td=" << degs.time;
    return ss.str();
  }

  int get_space_deg(const Point &p) const { return cell_deg.find(p)->second.space; }

  int get_space_deg(const Cell &c) const { return get_space_deg(c()); }

  int get_time_deg(const Point &p) const { return cell_deg.find(p)->second.time; }

  int get_time_deg(const Cell &c) const { return get_time_deg(c()); }

  void set_space_deg(const Point &p, int deg_space) {
    DegreePair old = cell_deg.find(p())->second;
    cell_deg[p] = DegreePair{short(deg_space), old.time};
    if (deg_space > max_deg.space) max_deg.space = deg_space;
  }

  void set_space_deg(const Cell &c, int deg) { set_space_deg(c(), deg); }

  void set_time_deg(const Point &p, int deg_time) {
    DegreePair old = cell_deg.find(p())->second;
    cell_deg[p] = DegreePair{old.space, short(deg_time)};
    if (deg_time > max_deg.time) max_deg.time = deg_time;
  }

  void set_time_deg(const Cell &c, int deg) { set_time_deg(c(), deg); }

  void increase_space_deg(const Cell &c) {
    DegreePair old = cell_deg.find(c())->second;
    cell_deg[c()] = DegreePair{short(old.space + 1), old.time};
    if (old.space + 1 > max_deg.space) max_deg.space = old.space + 1;
  }

  void decrease_space_deg(const Cell &c) {
    DegreePair old = cell_deg.find(c())->second;
    if (old.space > LOWEST_POSSIBLE_TIME_DEG) {
      cell_deg[c()] = DegreePair{short(old.space - 1), old.time};
    }
  }

  void increase_time_deg(const Cell &c) {
    DegreePair old = cell_deg.find(c())->second;
    cell_deg[c()] = DegreePair{old.space, short(old.time + 1)};
    if (old.time + 1 > max_deg.time) max_deg.time = old.time + 1;
  }

  void decrease_time_deg(const Cell &c) {
    DegreePair old = cell_deg.find(c())->second;
    if (old.time > LOWEST_POSSIBLE_TIME_DEG) {
      cell_deg[c()] = DegreePair{old.space, short(old.time - 1)};
    }
  }

  SpaceTimeMatrixGraph(const Mesh &stMesh, std::unique_ptr<IDoF> _dofs,
                       DegreePair defaultDegreePair, bool isDgInTime, bool single, LevelPair level,
                       const std::unordered_map<Point, DegreePair> &polyDist = {},
                       bool blockDiagonal = false);

  void update() override;

  bool IsSpaceTime() const override { return true; }

  ~SpaceTimeMatrixGraph() override = default;
};

#endif // SPACETIMEMATRIXGRAPH_HPP
