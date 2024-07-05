#ifndef ARGYRISDOFS_HPP
#define ARGYRISDOFS_HPP

#include "IDoF.hpp"

class ArgyrisDoF : public IDoF {
  short m;
public:
  explicit ArgyrisDoF(int M = 1, int n = 0, bool bnd = false) : IDoF(n, bnd), m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  std::vector<short> AccumulatedAllocationSizes(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const override { return "ArgyrisDoF"; }
};

#endif // ARGYRISDOFS_HPP
