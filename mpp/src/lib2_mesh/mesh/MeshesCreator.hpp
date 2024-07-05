#ifndef MESHESCREATOR_HPP
#define MESHESCREATOR_HPP

#include "Distribution.hpp"
#include "MeshSettings.hpp"
#include "Meshes.hpp"

class MeshesCreator {
private:
  // If no CoarseGeometry is set, 'meshName' is used to set up a CoarseGeometry
  //  If meshName is still "" for Create(), no CoarseGeometry will be produced
  std::string meshName = "";

  std::vector<double> timesteps{};

  MeshSettings settings;
public:
  MeshesCreator();

  explicit MeshesCreator(std::string meshName);

  explicit MeshesCreator(const MeshSettings &settings);

  const MeshSettings &Settings() const { return settings; }

  std::string GetMeshName() const { return meshName; }

  MeshesCreator WithMeshName(const std::string &meshName);

  MeshesCreator WithDirichletAt(std::vector<int> faceIndex);

  MeshesCreator WithNeumannAt(std::vector<int> faceIndex);

  MeshesCreator WithTimeRefinement(int timeRefinement);

  MeshesCreator WithFinalBarycentricRefinement();

  MeshesCreator WithShiftedCenter();

  MeshesCreator WithBalanceWeights(Distribution::Weights &&weights);

  MeshesCreator WithProcIds(std::vector<int> procIds);

  MeshesCreator WithLevel(int level);

  MeshesCreator WithPLevel(int pLevel);

  MeshesCreator WithCLevel(int cLevel);

  MeshesCreator WithDistribute(const std::string &distName);

  MeshesCreator WithoutDistribute();

  MeshesCreator WithCoarseGeometry(std::shared_ptr<CoarseGeometry> cGeo);

  MeshesCreator WithTimesteps(const std::vector<double> &ts);

  MeshesCreator WithDistribution(std::shared_ptr<Distribution> dist);

  MeshesCreator WithIncludeVector(std::vector<int> include);

  MeshesCreator WithSubdomainSetter(std::function<short(const Cell &)> sdSetter);

#ifdef USE_DATAMESH
  MeshesCreator WithCellData(DataContainer data);

  MeshesCreator WithVertexData(DataContainer data);
#endif

  Meshes *Create();

  std::unique_ptr<Meshes> CreateUnique();

  std::shared_ptr<Meshes> CreateShared();
};

#endif // MESHESCREATOR_HPP
