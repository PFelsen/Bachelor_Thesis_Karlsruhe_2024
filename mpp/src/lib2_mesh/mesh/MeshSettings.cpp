#include "MeshSettings.hpp"

#include "Config.hpp"

void MeshSettings::readConfig() {
  Config::Get("level", fineLevel);
  distributeLevel = fineLevel;
  Config::Get("plevel", distributeLevel);
  coarseLevel = distributeLevel;
  Config::Get("clevel", coarseLevel);
  Config::Get("Distribution", distributionName);
  Config::Get("Scales", timeRefinement);
  Config::Get("MeshesVerbose", verbose);
}

void MeshSettings::checkLevels() {
  if (coarseLevel > fineLevel) {
    coarseLevel = fineLevel;
    Warning("CoarseLevel set to fineLevel= " + std::to_string(fineLevel))
  }
  if (distributeLevel > fineLevel) {
    distributeLevel = fineLevel;
    Warning("DistributeLevel set to fineLevel= " + std::to_string(fineLevel))
  }
  coarseLevel = std::max(coarseLevel, distributeLevel);
}

void MeshSettings::initCoarseGeometry(const std::string &meshName) {
  if (!coarseGeometry) {
    if (meshName.empty()) { THROW("No meshName set.") }
    coarseGeometry =
        CreateCoarseGeometryShared(meshName, VertexDataList(vertexData), CellDataList(cellData));
  }
}
