#ifndef RTDOF_HPP
#define RTDOF_HPP

#include "IDoF.hpp"
#include "RTNodalPoints.hpp"

// TODO: Extend RTNodalPoints and RTDoF to all cell types

class RTDoF : public IDoF {
  int m;
  int order;
public:
  explicit RTDoF(int order, int M = 1, int n = 0, bool bnd = false) :
      IDoF(n, bnd), m(M), order(order) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const override { return "RaviartThomasDoF (order=" + std::to_string(order) + ")"; }
};

#endif // RTDOF_HPP
