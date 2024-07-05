#ifndef MESHES_HPP
#define MESHES_HPP

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

#include "CoarseGeometry.hpp"
#include "Distribution.hpp"
#include "LevelPair.hpp"
#include "Mesh.hpp"
#include "MeshSettings.hpp"
#include "Parallel.hpp"

class Meshes {
private:
  struct LevelPairTruncHash {
    std::size_t operator()(const LevelPair &levels) const {
      return levels.space * 10000 + levels.time + 100 * levels.commSplit;
    }
  };

  struct LevelPairTruncEqual {
    bool operator()(const LevelPair &l1, const LevelPair &l2) const {
      return (l1.space == l2.space) && (l1.time == l2.time) && (l1.commSplit == l2.commSplit);
    }
  };
public:
  using MeshesContainer =
      std::unordered_map<LevelPair, Mesh, LevelPairTruncHash, LevelPairTruncEqual>;

  /**
   * @brief Retrieve a mesh
   *
   * @param lp the level. The caller has to enforce following restrictions:
   *  - space must be >= settings.distributeLevel
   *  - lp.time must be -1 for space only meshes
   *  - lp.time must be >= settings.distributeLevel for ST meshes.
   * @return const Mesh&
   */
  const Mesh &operator[](const LevelPair &lp) const;

  /**
   * @brief Retrieve or create a new Mesh
   *
   * @param level the level.
   *  - Must be >= settings.distributeLevel and <= settings.fineLevel
   * @return const Mesh&
   */
  const Mesh &operator[](int level) const {
    if (isSTMeshes) { return (*this)[LevelPair{level, level + settings.timeRefinement, 0, 0}]; }
    return (*this)[LevelPair{level, -1, 0, 0}];
  }

  /**
   * @brief Removes a level-mesh pair
   *
   * @param lp the level
   * @return true if the pair has been removed
   * @return false otherwise
   */
  bool Erase(const LevelPair &lp) { return meshes.erase(lp); }

  /**
   * @brief Returns the number of meshes in this container
   *
   * @return std::size_t the number of meshes
   */
  std::size_t Size() const { return meshes.size(); }

  /**
   * @brief Returns whether this container owns a mesh with the given level
   *
   * @return bool true, if this container owns the requested mesh. false
   * otherwise
   */
  bool Contains(const LevelPair &lp) const { return meshes.contains(lp); }

  /**
   * @brief Removes all meshes in this container
   */
  void Clear() { meshes.clear(); }

  const MeshSettings &Settings() const { return settings; }

  std::string Name() const { return settings.coarseGeometry->Name(); }

  const CoarseGeometry &Domain() const { return *settings.coarseGeometry; }

  /**
   * @brief Get the Distribution of meshes with the given commSplit
   *
   * @param commSplit
   * @return std::shared_ptr<Distribution>, might be a nullptr
   */
  std::shared_ptr<Distribution> GetDistribution(int commSplit = 0) const {
    return settings.distributions[commSplit];
  }

  int dim() const { return fine().dim(); }

  LevelPair FineLevel() const {
    if (isSTMeshes) {
      return LevelPair{settings.fineLevel, settings.fineLevel + settings.timeRefinement, 0, 0};
    }
    return LevelPair{settings.fineLevel, -1, 0, 0};
  }

  LevelPair CoarseLevel() const {
    if (isSTMeshes) {
      return LevelPair{settings.coarseLevel, settings.coarseLevel + settings.timeRefinement, 0, 0};
    }
    return LevelPair{settings.coarseLevel, -1, 0, 0};
  }

  LevelPair PLevel() const {
    if (isSTMeshes) {
      return LevelPair{settings.distributeLevel, settings.distributeLevel + settings.timeRefinement,
                       0, 0};
    }
    return LevelPair{settings.distributeLevel, -1, 0, 0};
  }

  int Level() const { return settings.fineLevel; }

  int pLevel() const { return settings.distributeLevel; }

  const Mesh &coarse() const { return (*this)[CoarseLevel()]; }

  const Mesh &fine() const { return (*this)[FineLevel()]; }

  LevelPair AdaptLevels(LevelPair levels) const {
    if (levels.space < 0 || levels.space > settings.fineLevel) {
      levels.space = settings.fineLevel;
    }

    if (!isSTMeshes) {
      levels.time = -1;
    } else if (levels.time < 0 || levels.time > settings.fineLevel) {
      levels.time = settings.fineLevel;
    }

    return levels;
  }

  void PrintInfo() const;
protected:
  mutable MeshSettings settings;
  mutable MeshesContainer meshes;

  const bool isSTMeshes;

  friend class MeshesCreator;

  explicit Meshes(MeshSettings &&settings);
};

#endif // MESHES_HPP
