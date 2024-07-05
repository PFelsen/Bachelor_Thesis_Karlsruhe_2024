#ifndef MESHSETTINGS_HPP
#define MESHSETTINGS_HPP

#include <functional>
#include <memory>
#include <unordered_map>

#include "CoarseGeometry.hpp"
#include "Distribution.hpp"

struct MeshSettings {
  friend class MeshesCreator;

  int fineLevel = 0;
  int distributeLevel = 0;
  int coarseLevel = 0;

  int timeRefinement = 0;

  bool barycentricRefinement = false;
  bool shiftedCenter = false;

  int verbose = 1;

  // Note: mesh name is stored in coarse geometry
  std::shared_ptr<CoarseGeometry> coarseGeometry = nullptr;

  std::function<short(const Cell &)> sdSetter = nullptr;

  std::string distributionName = "RCB";

  using CommSplitDistributionMap = std::unordered_map<int, std::shared_ptr<Distribution>>;
  CommSplitDistributionMap distributions;

  std::shared_ptr<Distribution::Weights> weights = std::make_shared<Distribution::Weights>();
  std::vector<int> include = {};

  VertexDataList vertexData{};
  CellDataList cellData{};
private:
  void readConfig();

  void checkLevels();

  void initCoarseGeometry(const std::string &meshName);
};

#endif // MESHSETTINGS_HPP
