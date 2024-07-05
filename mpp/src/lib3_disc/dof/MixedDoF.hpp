#ifndef MIXEDDOF_HPP
#define MIXEDDOF_HPP

#include "IDoF.hpp"

struct SPData {
  int index;
  int shift;
};

struct NPData {
  int index;
};

class MixedDoF : public IDoF {
protected:
  std::vector<std::unique_ptr<IDoF>> dofs;

  mutable std::map<CELLTYPE, int> numNodalPoints{};
  mutable std::map<CELLTYPE, vector<short>> dofSizes{};
  mutable std::map<CELLTYPE, vector<vector<int>>> idNodalPointsOnFace{};
  mutable std::map<CELLTYPE, vector<vector<int>>> idNodalPointsOnEdge{};

  mutable std::map<CELLTYPE, int> numStoragePoints{};
  mutable std::map<CELLTYPE, vector<vector<int>>> idStoragePointsOnFace{};
  mutable std::map<CELLTYPE, vector<vector<int>>> idStoragePointsOnEdge{};
  mutable std::map<CELLTYPE, vector<short>> allocationsSizes{};
  mutable std::map<CELLTYPE, vector<short>> accumulatedAllocationSizes{};
  /**
   * storagePointData[n][i] contains the index and shift information of the i-th shape
   * function of the n-th shape with respect to ALL rows (not only the rows needed for
   * the n-th shape) on the Cell
   */
  mutable std::map<CELLTYPE, vector<vector<SPData>>> storagePointData{};
  mutable std::map<CELLTYPE, vector<vector<NPData>>> nodalPointData{};

  void initNodalPointData(const Cell &c) const;

  void initNodalPointOnFaceData(const Cell &c) const;

  void initNodalPointOnEdgeData(const Cell &c) const;

  void initStoragePointData(const Cell &c) const;

  void initStoragePointOnFaceData(const Cell &c) const;

  void initStoragePointOnEdgeData(const Cell &c) const;
public:
  static const MixedDoF &Cast(const IDoF &dof);

  MixedDoF(std::vector<std::unique_ptr<IDoF>> &&_dofs, int n = 0, bool b = false);

  short NumberOfComponents(int n) const;

  const vector<vector<SPData>> &StoragePointData(const Cell &c) const;

  const vector<SPData> &StoragePointData(int n, const Cell &c) const;

  const SPData &StoragePointData(int n, const Cell &c, int i) const;

  short NumberOfDoFs() const override { return dofs.size(); }

  short NumberOfNodalPoints(int n, const Cell &c) const;

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(int n, const Cell &c) const;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(int n, const Cell &c) const;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(int n, const Cell &c, int faceId) const;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(int n, const Cell &c, int faceId, int k) const;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  std::vector<Point> GetNodalPointsOnFace(int n, const Cell &c, int faceId) const;

  short NumberOfNodalPointsOnEdge(int n, const Cell &c, int edgeId) const;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(int n, const Cell &c, int edgeId, int k) const;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  std::vector<Point> GetNodalPointsOnEdge(int n, const Cell &c, int edgeId) const;

  short NumberOfStoragePoints(int n, const Cell &c) const;

  virtual short NumberOfStoragePoints(const Cell &c) const override;

  std::vector<Point> GetStoragePoints(int n, const Cell &c) const;

  virtual std::vector<Point> GetStoragePoints(const Cell &c) const override;

  std::vector<short> AllocationSizesAtStoragePoints(int n, const Cell &c) const;

  virtual std::vector<short> AllocationSizesAtStoragePoints(const Cell &c) const override;

  std::vector<short> AccumulatedAllocationSizes(int n, const Cell &c) const;

  short NumberOfStoragePointsOnFace(int n, const Cell &c, int faceId) const;

  virtual short NumberOfStoragePointsOnFace(const Cell &c, int faceId) const override;

  short IdOfStoragePointOnFace(int n, const Cell &c, int faceId, int k) const;

  virtual short IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const override;

  std::vector<Point> GetStoragePointsOnFace(int n, const Cell &c, int faceId) const;

  short NumberOfStoragePointsOnEdge(int n, const Cell &c, int edgeId) const;

  virtual short NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfStoragePointOnEdge(int n, const Cell &c, int edgeId, int k) const;

  virtual short IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const override;

  std::vector<Point> GetStoragePointsOnEdge(int n, const Cell &c, int edgeId) const;

  virtual std::string Name() const override;
};

#endif // MIXEDDOF_HPP
