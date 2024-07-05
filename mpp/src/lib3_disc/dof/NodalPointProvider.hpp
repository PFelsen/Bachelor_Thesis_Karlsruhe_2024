#ifndef SPACETIME_NODALPOINTPROVIDER_HPP
#define SPACETIME_NODALPOINTPROVIDER_HPP

#include "Celltype.hpp"
#include "GaussLobattoNodalPoints.hpp"
#include "Point.hpp"

#include "Cell.hpp"

class NodalPointProvider {
public:
  virtual std::vector<Point> GetNodalPoints(CELLTYPE ctype, int degree) const = 0;
  virtual ~NodalPointProvider() = default;
};

class EquidistantNodalPointProvider : public NodalPointProvider {
public:
  std::vector<Point> GetNodalPoints(CELLTYPE ctype, int degree) const override;
  ~EquidistantNodalPointProvider() override = default;
};

class GaussLobattoNodalPointProvider : public NodalPointProvider {
public:
  std::vector<Point> GetNodalPoints(CELLTYPE ctype, int degree) const override;
  ~GaussLobattoNodalPointProvider() override = default;
};

int EquidistantNodalPointCount(CELLTYPE type, int degree);

int EquidistantNodalPointCountOnFace(CELLTYPE type, int degree, int face);

int EquidistantNodalPointIdOnFace(CELLTYPE type, int degree, int face, int k);

int EquidistantNodalPointCountOnEdge(CELLTYPE type, int degree, int edge);

int EquidistantNodalPointIdOnEdge(CELLTYPE type, int degree, int edge, int k);

void EquidistantNodalPoints(const Cell &c, int degree, std::vector<Point> &z);

std::vector<Point> EquidistantNodalPoints(const Cell &c, int degree);

#endif // SPACETIME_NODALPOINTPROVIDER_HPP
