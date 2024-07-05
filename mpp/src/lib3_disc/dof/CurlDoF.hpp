#ifndef CURLDOF_HPP
#define CURLDOF_HPP

#include "IDoF.hpp"

class CurlDoF : public IDoF {
  short m;
  int degree;
public:
  explicit CurlDoF(int degree, short M = 1) : m(M), degree(degree) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const override { return "CurlDoF (degree=" + std::to_string(degree) + ")"; }
};

// The following DoFs need to be implemented in CurlDoF above

class Curl2TetDoF : public IDoF {
  short m;
public:
  Curl2TetDoF(short M = 1) : m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const { return "Curl2TetDoF"; }

  int SmallestFaceCorner(const Cell &c, int i) const;
};

class Curl2DoF : public IDoF {
  short m;

  // remapping of face corners local numbering, need to coincide the DOFs
  Point FaceCor(const Cell &c, int i, int j) const;
public:
  Curl2DoF(short M = 1) : m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const { return "Curl2DoF"; }
};

#endif // CURLDOF_HPP
