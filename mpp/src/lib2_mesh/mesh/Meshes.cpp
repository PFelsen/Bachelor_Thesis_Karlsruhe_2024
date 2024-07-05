#include "Meshes.hpp"

#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Distribution.hpp"
#include "LevelPair.hpp"
#include "Mesh.hpp"
#include "MeshSettings.hpp"
#include "Parallel.hpp"

#ifndef USE_SPACETIME
#define USE_SPACETIME 0
#endif

namespace {

Distribution &getDistribution(MeshSettings &settings, int commSplit) {
  auto &distPtr = settings.distributions[commSplit];
  if (distPtr == nullptr) {
    distPtr =
        std::make_shared<Distribution>(settings.distributionName, settings.weights, commSplit);
  }
  return *distPtr;
}

Mesh &refineAndSave(Mesh &baseMesh, Meshes::MeshesContainer &container, int toSpaceLevel,
                    int toTimeLevel = -1) {
  std::reference_wrapper<Mesh> currentMesh = baseMesh;

  std::vector<std::reference_wrapper<Mesh>> refineInSpaceCandidates;
  const int timeRefinement = toTimeLevel - baseMesh.Level().time;
  refineInSpaceCandidates.reserve(timeRefinement + 1);
  refineInSpaceCandidates.push_back(currentMesh);

  // Refine in time
  for (LevelPair lp = baseMesh.Level().NextInTime(); lp.time <= toTimeLevel; ++lp.time) {
    const auto &[iter, _] = container.emplace(std::piecewise_construct, std::forward_as_tuple(lp),
                                              std::forward_as_tuple(currentMesh.get(), true, true));
    Mesh &refinedMesh = iter->second;
    refinedMesh.InitTValues();
    currentMesh = refinedMesh;
    refineInSpaceCandidates.push_back(refinedMesh);
  }

  // Refine in space
  for (const auto &mesh : refineInSpaceCandidates) {
    currentMesh = mesh;
    for (LevelPair lp = mesh.get().Level().NextInSpace(); lp.space <= toSpaceLevel; ++lp.space) {
      const auto &[iter, _] = container.emplace(std::piecewise_construct, std::forward_as_tuple(lp),
                                                std::forward_as_tuple(currentMesh, true));
      Mesh &refinedMesh = iter->second;
      if (toTimeLevel >= 0) { refinedMesh.InitTValues(); }
      currentMesh = refinedMesh;
    }
  }

  return currentMesh;
}

Mesh &createUndistributedMappedMesh(const MeshSettings &settings, int commSplit,
                                    Meshes::MeshesContainer &container) {
  if (settings.distributeLevel == 0) {
    // Emplace (mapped) lowest level mesh
    // TODO: find an other way to avoid invoking the copy constructor
    if (settings.coarseGeometry->HasMapGeometry()) {
      // TODO: The last "true" in the Mesh constructor stands for mapGeometry.
      // Currently ONLY THE MESH ON THE DISTRIBUTED level is mapped,
      // which introduces larger boundary errors on fine levels,
      // but is used this way in DGWAVE
      return container
          .emplace(std::piecewise_construct, std::forward_as_tuple(LevelPair{0, -1, 0, commSplit}),
                   std::forward_as_tuple(Mesh(settings, commSplit), false, false, true))
          .first->second;
    }
    return container
        .emplace(std::piecewise_construct, std::forward_as_tuple(LevelPair{0, -1, 0, commSplit}),
                 std::forward_as_tuple(settings, commSplit))
        .first->second;
  }

  // Discard < distributeLevel meshes
  std::unique_ptr<Mesh> tmpMesh = std::make_unique<Mesh>(settings, commSplit);
  for (int nextSpaceLevel = 1; nextSpaceLevel < settings.distributeLevel; ++nextSpaceLevel) {
    tmpMesh = std::make_unique<Mesh>(*tmpMesh, true);
  }

  // Emplace mesh at distribute level
  return container
      .emplace(std::piecewise_construct, std::forward_as_tuple(tmpMesh->Level().NextInSpace()),
               std::forward_as_tuple(*tmpMesh, true, false,
                                     settings.coarseGeometry->HasMapGeometry()))
      .first->second;
}

void createBaseMeshes(Meshes::MeshesContainer &meshes, MeshSettings &settings, bool isSTMeshes,
                      int commSplit) {
  const int spaceMeshCount = settings.fineLevel - settings.distributeLevel + 1;
  const int spaceTimeMeshCount =
      (settings.fineLevel - settings.distributeLevel + 1)
      * (settings.fineLevel - settings.distributeLevel + 1 + settings.timeRefinement);

  if (isSTMeshes) {
    meshes.reserve(spaceMeshCount + spaceTimeMeshCount);
  } else {
    meshes.reserve(spaceMeshCount);
  }

  // Create (distributed) space mesh
  std::reference_wrapper<Mesh> distributeLevelSpaceMesh =
      createUndistributedMappedMesh(settings, commSplit, meshes);
  if (!isSTMeshes && PPM->IsParallel(commSplit)) {
    getDistribution(settings, commSplit).DistributeMesh(distributeLevelSpaceMesh);
  }
  // Create meshes up to fine mesh
  refineAndSave(distributeLevelSpaceMesh, meshes, settings.fineLevel);

  if (isSTMeshes) {
    std::reference_wrapper<Mesh> currentMesh = distributeLevelSpaceMesh;
    // Discard ST Meshes with timeLevel < distributionLevel
    std::unique_ptr<Mesh> tmpMesh;
    for (int nextTimeLevel = 0; nextTimeLevel < settings.distributeLevel; ++nextTimeLevel) {
      tmpMesh = std::make_unique<Mesh>(currentMesh, true, true);
      currentMesh = *tmpMesh;
    }

    // Create distributed ST mesh
    currentMesh = meshes
                      .emplace(std::piecewise_construct,
                               std::forward_as_tuple(currentMesh.get().Level().NextInTime()),
                               std::forward_as_tuple(currentMesh, true, true))
                      .first->second;
    if (PPM->IsParallel(commSplit)) {
      getDistribution(settings, commSplit).DistributeMesh(currentMesh);
    }

    // Create ST meshes up to fine mesh
    refineAndSave(currentMesh, meshes, settings.fineLevel,
                  settings.fineLevel + settings.timeRefinement);
  }
}

} // namespace

Meshes::Meshes(MeshSettings &&settings) :
    settings(settings), isSTMeshes(!settings.coarseGeometry->timeSteps.empty() && USE_SPACETIME) {
  createBaseMeshes(meshes, this->settings, isSTMeshes, 0);
}

const Mesh &Meshes::operator[](const LevelPair &lp) const {
  auto it = meshes.find(lp);
  if (it != std::end(meshes)) { return it->second; }

  // Find closest to lp
  LevelPair found{-1, -1, 0, lp.commSplit};
  int distance = std::numeric_limits<int>::max();
  for (const auto &[level, _] : meshes) {
    if (lp.commSplit != level.commSplit) { continue; }

    const int distanceToLp = lp.space - level.space + lp.time - level.time;

    if (level.space <= lp.space && level.time <= lp.time && distanceToLp < distance) {
      found = level;
      distance = distanceToLp;
    }
  }

  if (found.space == -1) {
    // Create base meshes for new commsplit if no base has been found
    createBaseMeshes(meshes, settings, isSTMeshes, lp.commSplit);
    // Try again after base meshes were created
    return (*this)[lp];
  }

  // Create requested mesh from found base
  if (isSTMeshes) {
    meshes.reserve(meshes.size() + (lp.space - found.space + 1) * (lp.time - found.time + 1));
    return refineAndSave(meshes.at(found), meshes, lp.space, lp.time);
  }

  meshes.reserve(meshes.size() + lp.space - found.space);
  std::reference_wrapper<Mesh> base = meshes.at(found);
  if (found.space < settings.distributeLevel && lp.space >= settings.distributeLevel) {
    // Map and distribute if needed
    const Mesh &preDistributedMesh = refineAndSave(base, meshes, settings.distributeLevel - 1);
    base = meshes
               .emplace(std::piecewise_construct,
                        std::forward_as_tuple(preDistributedMesh.Level().NextInSpace()),
                        std::forward_as_tuple(preDistributedMesh, true, false,
                                              settings.coarseGeometry->HasMapGeometry()))
               .first->second;
    if (PPM->IsParallel(lp.commSplit)) {
      getDistribution(settings, lp.commSplit).DistributeMesh(base);
    }
  }
  return refineAndSave(base, meshes, lp.space);
}

void Meshes::PrintInfo() const {
  const int minLevel = settings.fineLevel - settings.verbose;
  for (const auto &[level, mesh] : meshes) {
    if (level.space >= minLevel) {
      if (!isSTMeshes || level.time >= minLevel + settings.timeRefinement) { mesh.PrintInfo(); }
    }
  }
}
