#ifndef DGDOF_HPP
#define DGDOF_HPP

#include "IDoF.hpp"
#include "NodalPointProvider.hpp"

class DGDoF : public IDoF {
  int cellStandardDegree;
  int m;
  std::unordered_map<Point, vector<short>> cellDegree;
public:
  explicit DGDoF(int cellStandardDegree, int M) :
      IDoF(0, false), m(M), cellStandardDegree(cellStandardDegree) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  short NumberOfStoragePoints(const Cell &c) const override;

  std::vector<Point> GetStoragePoints(const Cell &c) const override;

  std::vector<short> AllocationSizesAtStoragePoints(const Cell &c) const override;

  short NumberOfStoragePointsOnFace(const Cell &c, int faceId) const override;

  short IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const override { return "DGDoF"; }

  int TypeDoF(int) const override { return CELL; }

  int get_cell_deg(const Cell &) const override { return cellStandardDegree; }

  void set_cell_deg(const Cell &c, int deg) override;
};

#endif // DGDOF_HPP
