#ifndef SERENDIPITYDOF_HPP
#define SERENDIPITYDOF_HPP

#include <map>
#include "IDoF.hpp"
#include "NodalPointProvider.hpp"

class SerendipityDoF : public IDoF {
protected:
  short m;
public:
  explicit SerendipityDoF(short M = 1, int n = 0, bool bnd = false) : IDoF(n, bnd), m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const override { return "SerendipityDoF"; }
};

#endif // SERENDIPITYDOF_HPP
