#include "VectorMatrixBase.hpp"

template<>
const IDiscretization &VectorMatrixBase::GetDiscT<double, SpaceDimension, TimeDimension>() const {
  return *disc;
}

#ifdef BUILD_IA

template<>
const IAIDiscretization &
VectorMatrixBase::GetDiscT<IAInterval, SpaceDimension, TimeDimension>() const {
  if (iadisc) return *iadisc;
  THROW("IADiscretization not set!")
}

#endif

std::string VectorMatrixBase::DoFsName() const { return GetDoF().Name(); }

short VectorMatrixBase::NumberOfDoFs() const { return GetDoF().NumberOfDoFs(); }

short VectorMatrixBase::NumberOfComponents() const { return GetDoF().NumberOfComponents(); }

short VectorMatrixBase::NumberOfNodalPoints(const Cell &c) const {
  return GetDoF().NumberOfNodalPoints(c);
}

short VectorMatrixBase::NumberOfNodalPoints(int n, const Cell &c) const {
  return MixedDoF::Cast(GetDoF()).NumberOfNodalPoints(n, c);
}

short VectorMatrixBase::NumberOfNodalPoints(const int sDeg, const int tDeg) const {
  return GetDoF().NumberOfNodalPoints(sDeg, tDeg);
}

std::vector<Point> VectorMatrixBase::GetNodalPoints(const Cell &c) const {
  return GetDoF().GetNodalPoints(c);
}

std::vector<Point> VectorMatrixBase::GetNodalPoints(int n, const Cell &c) const {
  return MixedDoF::Cast(GetDoF()).GetNodalPoints(n, c);
}

std::vector<short> VectorMatrixBase::DoFSizesAtNodalPoints(const Cell &c) const {
  return GetDoF().DoFSizesAtNodalPoints(c);
}

std::vector<short> VectorMatrixBase::DoFSizesAtNodalPoints(int n, const Cell &c) const {
  return MixedDoF::Cast(GetDoF()).DoFSizesAtNodalPoints(n, c);
}

short VectorMatrixBase::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const {
  return GetDoF().NumberOfNodalPointsOnFace(c, faceId);
}

short VectorMatrixBase::NumberOfNodalPointsOnFace(int n, const Cell &c, int faceId) const {
  return MixedDoF::Cast(GetDoF()).NumberOfNodalPointsOnFace(n, c, faceId);
}

short VectorMatrixBase::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  return GetDoF().IdOfNodalPointOnFace(c, faceId, k);
}

short VectorMatrixBase::IdOfNodalPointOnFace(int n, const Cell &c, int faceId, int k) const {
  return MixedDoF::Cast(GetDoF()).IdOfNodalPointOnFace(n, c, faceId, k);
}

std::vector<Point> VectorMatrixBase::GetNodalPointsOnFace(const Cell &c, int faceId) const {
  return GetDoF().GetNodalPointsOnFace(c, faceId);
}

std::vector<Point> VectorMatrixBase::GetNodalPointsOnFace(int n, const Cell &c, int faceId) const {
  return MixedDoF::Cast(GetDoF()).GetNodalPointsOnFace(n, c, faceId);
}

short VectorMatrixBase::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const {
  return GetDoF().NumberOfNodalPointsOnEdge(c, edgeId);
}

short VectorMatrixBase::NumberOfNodalPointsOnEdge(int n, const Cell &c, int edgeId) const {
  return MixedDoF::Cast(GetDoF()).NumberOfNodalPointsOnEdge(n, c, edgeId);
}

short VectorMatrixBase::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  return GetDoF().IdOfNodalPointOnEdge(c, edgeId, k);
}

short VectorMatrixBase::IdOfNodalPointOnEdge(int n, const Cell &c, int edgeId, int k) const {
  return MixedDoF::Cast(GetDoF()).IdOfNodalPointOnEdge(n, c, edgeId, k);
}

std::vector<Point> VectorMatrixBase::GetNodalPointsOnEdge(const Cell &c, int edgeId) const {
  return GetDoF().GetNodalPointsOnEdge(c, edgeId);
}

std::vector<Point> VectorMatrixBase::GetNodalPointsOnEdge(int n, const Cell &c, int edgeId) const {
  return MixedDoF::Cast(GetDoF()).GetNodalPointsOnEdge(n, c, edgeId);
}