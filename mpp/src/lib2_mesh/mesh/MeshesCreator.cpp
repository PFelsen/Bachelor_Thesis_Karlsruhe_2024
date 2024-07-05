#include "MeshesCreator.hpp"

#include <memory>
#include <utility>

#include "Distribution.hpp"

MeshesCreator::MeshesCreator() {
  Config::Get("Mesh", meshName);
  settings.readConfig();
}

MeshesCreator::MeshesCreator(std::string _meshName) : meshName(std::move(_meshName)) {
  settings.readConfig();
}

MeshesCreator::MeshesCreator(const MeshSettings &settings) :
    settings(settings), meshName(settings.coarseGeometry->Name()) {}

MeshesCreator MeshesCreator::WithMeshName(const std::string &_meshName) {
  meshName = _meshName;
  return *this;
}

MeshesCreator MeshesCreator::WithDirichletAt(std::vector<int> faceIndex) {
  Warning("MeshesCreator::WithDirichletAt(std::vector<int>) has not effect!") return *this;
}

MeshesCreator MeshesCreator::WithNeumannAt(std::vector<int> faceIndex) {
  Warning("MeshesCreator::WithNeumannAt(std::vector<int>) has not effect!") return *this;
}

MeshesCreator MeshesCreator::WithTimeRefinement(int timeRefinement) {
  settings.timeRefinement = timeRefinement;
  return *this;
}

MeshesCreator MeshesCreator::WithFinalBarycentricRefinement() {
  settings.barycentricRefinement = true;
  return *this;
}

MeshesCreator MeshesCreator::WithShiftedCenter() {
  settings.shiftedCenter = true;
  return *this;
}

MeshesCreator MeshesCreator::WithBalanceWeights(Distribution::Weights &&weights) {
  *settings.weights = std::move(weights);
  return *this;
}

MeshesCreator MeshesCreator::WithProcIds(std::vector<int> procIds) {
  Warning("MeshesCreator::WithProcIds(std::vector<int>) has not effect!") return *this;
}

MeshesCreator MeshesCreator::WithLevel(int level) {
  settings.fineLevel = level;
  return *this;
}

MeshesCreator MeshesCreator::WithPLevel(int pLevel) {
  settings.distributeLevel = pLevel;
  return *this;
}

MeshesCreator MeshesCreator::WithCLevel(int cLevel) {
  settings.coarseLevel = cLevel;
  return *this;
}

MeshesCreator MeshesCreator::WithDistribute(const std::string &distName) {
  settings.distributionName = distName;
  return *this;
}

MeshesCreator MeshesCreator::WithoutDistribute() {
  settings.distributionName = "NoDistribution";
  return *this;
}

MeshesCreator MeshesCreator::WithCoarseGeometry(std::shared_ptr<CoarseGeometry> cGeo) {
  settings.coarseGeometry = std::move(cGeo);
  return *this;
}

MeshesCreator MeshesCreator::WithTimesteps(const std::vector<double> &ts) {
  timesteps = ts;
  return *this;
}

MeshesCreator MeshesCreator::WithDistribution(std::shared_ptr<Distribution> dist) {
  settings.distributionName = dist->Name();
  settings.distributions[dist->CommSplit()] = std::move(dist);
  return *this;
}

MeshesCreator MeshesCreator::WithIncludeVector(std::vector<int> include) {
  settings.include = std::move(include);
  return *this;
}

MeshesCreator MeshesCreator::WithSubdomainSetter(std::function<short(const Cell &)> _sdSetter) {
  settings.sdSetter = std::move(_sdSetter);
  return *this;
}

#ifdef USE_DATAMESH
MeshesCreator MeshesCreator::WithCellData(DataContainer data) {
  settings.cellData = CellDataList(1, data);
  return *this;
}

MeshesCreator MeshesCreator::WithVertexData(DataContainer data) {
  settings.vertexData = VertexDataList(1, data);
  return *this;
}
#endif

Meshes *MeshesCreator::Create() {
  settings.checkLevels();
  settings.initCoarseGeometry(meshName);
  if (!timesteps.empty()) { settings.coarseGeometry->timeSteps = timesteps; }
  return new Meshes(std::move(settings));
}

std::unique_ptr<Meshes> MeshesCreator::CreateUnique() { return std::unique_ptr<Meshes>(Create()); }

std::shared_ptr<Meshes> MeshesCreator::CreateShared() { return std::shared_ptr<Meshes>(Create()); }
