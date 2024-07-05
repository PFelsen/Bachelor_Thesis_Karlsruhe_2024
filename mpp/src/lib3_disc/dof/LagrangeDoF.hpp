#ifndef LAGRANGEDOF_HPP
#define LAGRANGEDOF_HPP

#include <map>
#include "IDoF.hpp"
#include "NodalPointProvider.hpp"

class LagrangeDoF : public IDoF {
protected:
  short m;
  int degree;
public:
  explicit LagrangeDoF(int degree, int M = 1, int n = 0, bool bnd = false) :
      IDoF(n, bnd), m(M), degree(degree) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override {
    return EquidistantNodalPointCount(c.ReferenceType(), degree);
  }

  std::vector<Point> GetNodalPoints(const Cell &c) const override {
    return EquidistantNodalPoints(c, degree);
  }

  virtual std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>(NumberOfNodalPoints(c), m);
  }

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override {
    return EquidistantNodalPointCountOnFace(c.ReferenceType(), degree, faceId);
  }

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override {
    return EquidistantNodalPointIdOnFace(c.ReferenceType(), degree, faceId, k);
  }

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override {
    return EquidistantNodalPointCountOnEdge(c.ReferenceType(), degree, edgeId);
  }

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override {
    return EquidistantNodalPointIdOnEdge(c.ReferenceType(), degree, edgeId, k);
  }

  string Name() const override { return "LagrangeDoF (degree=" + std::to_string(degree) + ")"; }

  int Degree() const { return degree; }
};

#endif // LAGRANGEDOF_HPP
