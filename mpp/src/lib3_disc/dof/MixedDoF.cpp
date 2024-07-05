#include "MixedDoF.hpp"

void MixedDoF::initNodalPointData(const Cell &c) const {
  vector<vector<NPData>> data(dofs.size());
  vector<Point> nodalPoints = GetNodalPoints(c);
  vector<short> _dofSizes(nodalPoints.size(), 0);

  for (int n = 0; n < dofs.size(); ++n) {
    vector<Point> nodalPointsN = dofs[n]->GetNodalPoints(c);
    vector<short> allocSizesN = dofs[n]->DoFSizesAtNodalPoints(c);
    data[n].resize(nodalPointsN.size());
    for (int i = 0; i < nodalPointsN.size(); ++i) {
      int idx = std::distance(nodalPoints.begin(),
                              std::find(nodalPoints.begin(), nodalPoints.end(), nodalPointsN[i]));
      data[n][i].index = idx;
      _dofSizes[idx] += allocSizesN[i];
    }
  }
  nodalPointData[c.Type()] = data;
  dofSizes[c.Type()] = _dofSizes;
}

void MixedDoF::initNodalPointOnFaceData(const Cell &c) const {
  auto data_c = nodalPointData.find(c.Type());
  if (data_c == nodalPointData.end()) {
    initNodalPointData(c);
    data_c = nodalPointData.find(c.Type());
  }

  vector<vector<int>> w(c.Faces());
  for (int face = 0; face < c.Faces(); ++face) {
    w[face].resize(dofs[0]->NumberOfNodalPointsOnFace(c, face));
    for (int k = 0; k < w[face].size(); ++k)
      w[face][k] = dofs[0]->IdOfNodalPointOnFace(c, face, k);
    for (int n = 1; n < dofs.size(); ++n) {
      for (int k = 0; k < dofs[n]->NumberOfNodalPointsOnFace(c, face); ++k) {
        int idx = data_c->second[n][dofs[n]->IdOfNodalPointOnFace(c, face, k)].index;
        if (std::find(w[face].begin(), w[face].end(), idx) == w[face].end()) w[face].push_back(idx);
      }
    }
  }
  idNodalPointsOnFace[c.Type()] = w;
}

void MixedDoF::initNodalPointOnEdgeData(const Cell &c) const {
  auto data_c = nodalPointData.find(c.Type());
  if (data_c == nodalPointData.end()) {
    initNodalPointData(c);
    data_c = nodalPointData.find(c.Type());
  }

  vector<vector<int>> w(c.Edges());
  for (int edge = 0; edge < c.Edges(); ++edge) {
    w[edge].resize(dofs[0]->NumberOfNodalPointsOnEdge(c, edge));
    for (int k = 0; k < w[edge].size(); ++k)
      w[edge][k] = dofs[0]->IdOfNodalPointOnEdge(c, edge, k);
    for (int n = 1; n < dofs.size(); ++n) {
      for (int k = 0; k < dofs[n]->NumberOfNodalPointsOnEdge(c, edge); ++k) {
        int idx = data_c->second[n][dofs[n]->IdOfNodalPointOnEdge(c, edge, k)].index;
        if (std::find(w[edge].begin(), w[edge].end(), idx) == w[edge].end()) w[edge].push_back(idx);
      }
    }
  }
  idNodalPointsOnEdge[c.Type()] = w;
}

void MixedDoF::initStoragePointData(const Cell &c) const {
  vector<vector<SPData>> data(dofs.size());
  vector<Point> storagePoints = GetStoragePoints(c);
  vector<short> allocSizes(storagePoints.size(), 0);

  for (int n = 0; n < dofs.size(); ++n) {
    vector<Point> nodalPointsN = dofs[n]->GetStoragePoints(c);
    vector<short> allocSizesN = dofs[n]->AllocationSizesAtStoragePoints(c);
    data[n].resize(nodalPointsN.size());
    for (int i = 0; i < nodalPointsN.size(); ++i) {
      int idx =
          std::distance(storagePoints.begin(),
                        std::find(storagePoints.begin(), storagePoints.end(), nodalPointsN[i]));
      data[n][i].index = idx;
      data[n][i].shift = allocSizes[idx];
      allocSizes[idx] += allocSizesN[i];
    }
  }
  storagePointData[c.Type()] = data;
  allocationsSizes[c.Type()] = allocSizes;

  for (int i = 1; i < allocSizes.size(); ++i) {
    allocSizes[i] += allocSizes[i - 1];
  }
  accumulatedAllocationSizes[c.Type()] = allocSizes;
}

void MixedDoF::initStoragePointOnFaceData(const Cell &c) const {
  auto data_c = storagePointData.find(c.Type());
  if (data_c == storagePointData.end()) {
    initStoragePointData(c);
    data_c = storagePointData.find(c.Type());
  }

  vector<vector<int>> w(c.Faces());
  for (int face = 0; face < c.Faces(); ++face) {
    w[face].resize(dofs[0]->NumberOfStoragePointsOnFace(c, face));
    for (int k = 0; k < w[face].size(); ++k)
      w[face][k] = dofs[0]->IdOfStoragePointOnFace(c, face, k);
    for (int n = 1; n < dofs.size(); ++n) {
      for (int k = 0; k < dofs[n]->NumberOfStoragePointsOnFace(c, face); ++k) {
        int idx = data_c->second[n][dofs[n]->IdOfStoragePointOnFace(c, face, k)].index;
        if (std::find(w[face].begin(), w[face].end(), idx) == w[face].end()) w[face].push_back(idx);
      }
    }
  }
  idStoragePointsOnFace[c.Type()] = w;
}

void MixedDoF::initStoragePointOnEdgeData(const Cell &c) const {
  auto data_c = storagePointData.find(c.Type());
  if (data_c == storagePointData.end()) {
    initStoragePointData(c);
    data_c = storagePointData.find(c.Type());
  }

  vector<vector<int>> w(c.Edges());
  for (int edge = 0; edge < c.Edges(); ++edge) {
    w[edge].resize(dofs[0]->NumberOfStoragePointsOnEdge(c, edge));
    for (int k = 0; k < w[edge].size(); ++k)
      w[edge][k] = dofs[0]->IdOfStoragePointOnEdge(c, edge, k);
    for (int n = 1; n < dofs.size(); ++n) {
      for (int k = 0; k < dofs[n]->NumberOfStoragePointsOnEdge(c, edge); ++k) {
        int idx = data_c->second[n][dofs[n]->IdOfStoragePointOnEdge(c, edge, k)].index;
        if (std::find(w[edge].begin(), w[edge].end(), idx) == w[edge].end()) w[edge].push_back(idx);
      }
    }
  }
  idStoragePointsOnEdge[c.Type()] = w;
}

const MixedDoF &MixedDoF::Cast(const IDoF &dof) {
  if (typeid(dof) != typeid(MixedDoF)) THROW("Cannot cast to MixedDoF")
  return dynamic_cast<const MixedDoF &>(dof);
}

MixedDoF::MixedDoF(std::vector<std::unique_ptr<IDoF>> &&_dofs, int n, bool b) :
    IDoF(n, b), dofs(std::move(_dofs)) {
  if (dofs.size() < 1) THROW("Size of MixedDoF to small")
}

short MixedDoF::NumberOfComponents(int n) const { return dofs[n]->NumberOfComponents(); }

const vector<vector<SPData>> &MixedDoF::StoragePointData(const Cell &c) const {
  auto data_c = storagePointData.find(c.Type());
  if (data_c == storagePointData.end()) {
    initStoragePointData(c);
    data_c = storagePointData.find(c.Type());
  }
  return data_c->second;
}

const vector<SPData> &MixedDoF::StoragePointData(int n, const Cell &c) const {
  return StoragePointData(c)[n];
}

const SPData &MixedDoF::StoragePointData(int n, const Cell &c, int i) const {
  return StoragePointData(n, c)[i];
}

short MixedDoF::NumberOfNodalPoints(int n, const Cell &c) const {
  return dofs[n]->NumberOfNodalPoints(c);
}

short MixedDoF::NumberOfNodalPoints(const Cell &c) const {
  auto numNodalPoints_c = numNodalPoints.find(c.Type());
  if (numNodalPoints_c != numNodalPoints.end()) return numNodalPoints_c->second;
  std::vector<Point> np = GetNodalPoints(c);
  numNodalPoints[c.Type()] = np.size();
  return np.size();
}

std::vector<Point> MixedDoF::GetNodalPoints(int n, const Cell &c) const {
  return dofs[n]->GetNodalPoints(c);
}

std::vector<Point> MixedDoF::GetNodalPoints(const Cell &c) const {
  std::vector np = dofs[0]->GetNodalPoints(c);
  for (int n = 1; n < dofs.size(); ++n) {
    std::vector<Point> w = dofs[n]->GetNodalPoints(c);
    for (Point &p : w) {
      if (std::find(np.begin(), np.end(), p) == np.end()) { np.push_back(p); }
    }
  }
  return np;
}

std::vector<short> MixedDoF::DoFSizesAtNodalPoints(int n, const Cell &c) const {
  return dofs[n]->DoFSizesAtNodalPoints(c);
}

std::vector<short> MixedDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  auto dofSizes_c = dofSizes.find(c.Type());
  if (dofSizes_c == dofSizes.end()) {
    initNodalPointData(c);
    dofSizes_c = dofSizes.find(c.Type());
  }
  return dofSizes_c->second;
}

short MixedDoF::NumberOfNodalPointsOnFace(int n, const Cell &c, int faceId) const {
  return dofs[n]->NumberOfNodalPointsOnFace(c, faceId);
}

short MixedDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const {
  auto npOnFace_c = idNodalPointsOnFace.find(c.Type());
  if (npOnFace_c == idNodalPointsOnFace.end()) {
    initNodalPointOnFaceData(c);
    npOnFace_c = idNodalPointsOnFace.find(c.Type());
  }
  return npOnFace_c->second[faceId].size();
}

short MixedDoF::IdOfNodalPointOnFace(int n, const Cell &c, int faceId, int k) const {
  return dofs[n]->IdOfNodalPointOnFace(c, faceId, k);
}

short MixedDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const {
  auto idNpOnFace_c = idNodalPointsOnFace.find(c.Type());
  if (idNpOnFace_c == idNodalPointsOnFace.end()) {
    initNodalPointOnFaceData(c);
    idNpOnFace_c = idNodalPointsOnFace.find(c.Type());
  }
  return idNpOnFace_c->second[faceId][k];
}

std::vector<Point> MixedDoF::GetNodalPointsOnFace(int n, const Cell &c, int faceId) const {
  return dofs[n]->GetNodalPointsOnFace(c, faceId);
}

short MixedDoF::NumberOfNodalPointsOnEdge(int n, const Cell &c, int edgeId) const {
  return dofs[n]->NumberOfNodalPointsOnEdge(c, edgeId);
}

short MixedDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const {
  auto npOnEdge_c = idNodalPointsOnEdge.find(c.Type());
  if (npOnEdge_c == idNodalPointsOnEdge.end()) {
    initNodalPointOnEdgeData(c);
    npOnEdge_c = idNodalPointsOnEdge.find(c.Type());
  }
  return npOnEdge_c->second[edgeId].size();
}

short MixedDoF::IdOfNodalPointOnEdge(int n, const Cell &c, int edgeId, int k) const {
  return dofs[n]->IdOfNodalPointOnEdge(c, edgeId, k);
}

short MixedDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const {
  auto idNpOnEdge_c = idNodalPointsOnEdge.find(c.Type());
  if (idNpOnEdge_c == idNodalPointsOnEdge.end()) {
    initNodalPointOnEdgeData(c);
    idNpOnEdge_c = idNodalPointsOnEdge.find(c.Type());
  }
  return idNpOnEdge_c->second[edgeId][k];
}

std::vector<Point> MixedDoF::GetNodalPointsOnEdge(int n, const Cell &c, int edgeId) const {
  return dofs[n]->GetNodalPointsOnEdge(c, edgeId);
}

short MixedDoF::NumberOfStoragePoints(int n, const Cell &c) const {
  return dofs[n]->NumberOfStoragePoints(c);
}

short MixedDoF::NumberOfStoragePoints(const Cell &c) const {
  auto numNodalPoints_c = numStoragePoints.find(c.Type());
  if (numNodalPoints_c != numStoragePoints.end()) return numNodalPoints_c->second;
  std::vector<Point> sp = GetStoragePoints(c);
  numStoragePoints[c.Type()] = sp.size();
  return sp.size();
}

std::vector<Point> MixedDoF::GetStoragePoints(int n, const Cell &c) const {
  return dofs[n]->GetStoragePoints(c);
}

std::vector<Point> MixedDoF::GetStoragePoints(const Cell &c) const {
  std::vector sp = dofs[0]->GetStoragePoints(c);
  for (int n = 1; n < dofs.size(); ++n) {
    std::vector<Point> w = dofs[n]->GetStoragePoints(c);
    for (Point &p : w) {
      if (std::find(sp.begin(), sp.end(), p) == sp.end()) { sp.push_back(p); }
    }
  }
  return sp;
}

std::vector<short> MixedDoF::AllocationSizesAtStoragePoints(int n, const Cell &c) const {
  return dofs[n]->AllocationSizesAtStoragePoints(c);
}

std::vector<short> MixedDoF::AllocationSizesAtStoragePoints(const Cell &c) const {
  auto allocSizes_c = allocationsSizes.find(c.Type());
  if (allocSizes_c == allocationsSizes.end()) {
    initStoragePointData(c);
    allocSizes_c = allocationsSizes.find(c.Type());
  }
  return allocSizes_c->second;
}

std::vector<short> MixedDoF::AccumulatedAllocationSizes(int n, const Cell &c) const {
  return dofs[n]->AccumulatedAllocationSizes(c);
}

short MixedDoF::NumberOfStoragePointsOnFace(int n, const Cell &c, int faceId) const {
  return dofs[n]->NumberOfStoragePointsOnFace(c, faceId);
}

short MixedDoF::NumberOfStoragePointsOnFace(const Cell &c, int faceId) const {
  auto spOnFace_c = idStoragePointsOnFace.find(c.Type());
  if (spOnFace_c == idStoragePointsOnFace.end()) {
    initStoragePointOnFaceData(c);
    spOnFace_c = idStoragePointsOnFace.find(c.Type());
  }
  return spOnFace_c->second[faceId].size();
}

short MixedDoF::IdOfStoragePointOnFace(int n, const Cell &c, int faceId, int k) const {
  return dofs[n]->IdOfStoragePointOnFace(c, faceId, k);
}

short MixedDoF::IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const {
  auto idSpOnFace_c = idStoragePointsOnFace.find(c.Type());
  if (idSpOnFace_c == idStoragePointsOnFace.end()) {
    initStoragePointOnFaceData(c);
    idSpOnFace_c = idStoragePointsOnFace.find(c.Type());
  }
  return idSpOnFace_c->second[faceId][k];
}

std::vector<Point> MixedDoF::GetStoragePointsOnFace(int n, const Cell &c, int faceId) const {
  return dofs[n]->GetStoragePointsOnFace(c, faceId);
}

short MixedDoF::NumberOfStoragePointsOnEdge(int n, const Cell &c, int edgeId) const {
  return dofs[n]->NumberOfStoragePointsOnEdge(c, edgeId);
}

short MixedDoF::NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const {
  auto spOnEdge_c = idStoragePointsOnEdge.find(c.Type());
  if (spOnEdge_c == idStoragePointsOnEdge.end()) {
    initStoragePointOnEdgeData(c);
    spOnEdge_c = idStoragePointsOnEdge.find(c.Type());
  }
  return spOnEdge_c->second[edgeId].size();
}

short MixedDoF::IdOfStoragePointOnEdge(int n, const Cell &c, int edgeId, int k) const {
  return dofs[n]->IdOfStoragePointOnEdge(c, edgeId, k);
}

short MixedDoF::IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const {
  auto idSpOnEdge_c = idStoragePointsOnEdge.find(c.Type());
  if (idSpOnEdge_c == idStoragePointsOnEdge.end()) {
    initStoragePointOnEdgeData(c);
    idSpOnEdge_c = idStoragePointsOnEdge.find(c.Type());
  }
  return idSpOnEdge_c->second[edgeId][k];
}

std::vector<Point> MixedDoF::GetStoragePointsOnEdge(int n, const Cell &c, int edgeId) const {
  return dofs[n]->GetStoragePointsOnEdge(c, edgeId);
}

std::string MixedDoF::Name() const {
  string name = dofs[0]->Name();
  for (int n = 1; n < dofs.size(); ++n)
    name.append("; ").append(dofs[n]->Name());
  return name;
}
