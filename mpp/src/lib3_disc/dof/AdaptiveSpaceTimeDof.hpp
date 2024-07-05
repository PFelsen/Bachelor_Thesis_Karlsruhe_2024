#ifndef ADAPTIVESPACETIMEDOF_HPP
#define ADAPTIVESPACETIMEDOF_HPP

#include "IDoF.hpp"
#include "SpaceTimeMatrixGraph.hpp"

class AdaptiveSpaceTimeDof : public IDoF {
protected:
  const int LOWEST_POSSIBLE_TIME_DEG = 0;

  SpaceTimeMatrixGraph *matrixGraph;
public:
  AdaptiveSpaceTimeDof(const int lowest_possible_time_deg) :
      LOWEST_POSSIBLE_TIME_DEG(lowest_possible_time_deg){};

  short NumberOfStoragePoints(const Cell &c) const override { return 1; }

  std::vector<Point> GetStoragePoints(const Cell &c) const override {
    return std::vector<Point>{c()};
  }

  void SetMatrixGraph(SpaceTimeMatrixGraph *mGraph) { matrixGraph = mGraph; }

  const int GetLOWEST_POSSIBLE_TIME_DEG() { return LOWEST_POSSIBLE_TIME_DEG; }

  short NumberOfStoragePointsOnFace(const Cell &c, int faceId) const override { return 0; }

  short IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const override {
    THROW("Should never reach this function")
  }

  short NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const override { return 0; }

  short IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const override {
    THROW("Should never reach this function")
  }

  void communicate() override { matrixGraph->communicate(); }

  void SetDegree(Point p, DegreePair pair) override { matrixGraph->SetDegree(p, pair); }

  DegreePair GetDegree(Point p) const override { return matrixGraph->GetDegree(p); }

  DegreePair GetDegree(const Cell &c) const override { return matrixGraph->GetDegree(c()); }

  DegreePair GetMaxDegree() const override { return matrixGraph->GetMaxDegree(); }

  string Name(const Cell &c) const override { return matrixGraph->Name(c); }

  int get_space_deg(const Point &p) const override { return matrixGraph->get_space_deg(p); }

  int get_space_deg(const Cell &c) const override { return matrixGraph->get_space_deg(c()); }

  int get_time_deg(const Point &p) const override { return matrixGraph->get_time_deg(p); }

  int get_time_deg(const Cell &c) const override { return matrixGraph->get_time_deg(c()); }

  void set_space_deg(const Point &p, int deg_space) override {
    matrixGraph->set_space_deg(p, deg_space);
  }

  void set_space_deg(const Cell &c, int deg) override { matrixGraph->set_space_deg(c(), deg); }

  void set_time_deg(const Point &p, int deg_time) override {
    matrixGraph->set_time_deg(p, deg_time);
  }

  void set_time_deg(const Cell &c, int deg) override { matrixGraph->set_time_deg(c, deg); }

  void increase_space_deg(const Cell &c) override { matrixGraph->increase_space_deg(c); }

  void decrease_space_deg(const Cell &c) override { matrixGraph->decrease_space_deg(c); }

  void increase_time_deg(const Cell &c) override { matrixGraph->increase_time_deg(c); }

  void decrease_time_deg(const Cell &c) override { matrixGraph->decrease_time_deg(c); }
};


#endif // SPACETIME_ADAPTIVESPACETIMEDOF_HPP
