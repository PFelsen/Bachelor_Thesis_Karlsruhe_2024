#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#include <memory>
#include "Config.hpp"
#include "Meshes.hpp"
#include "MeshesCreator.hpp"

enum class COMPONENT {
  V_X,
  V_Y,
  V_Z,
  P0,
  P1,
  P2,
  P3,
  P4,
  P5,
  S0_x,
  S0_y,
  S0_xy,
  S1_x,
  S1_y,
  S1_xy,
  S2_x,
  S2_y,
  S2_xy,
  S3_x,
  S3_y,
  S3_xy,
  S4_x,
  S4_y,
  S4_xy
};

std::string to_string(COMPONENT comp);

bool isVelocity(COMPONENT comp);

int toVelocityIndex(COMPONENT comp);

COMPONENT GetDampingComponent(int i);

class IProblem {
protected:
  int verbose = 1;

  std::string meshesName = "";

  std::shared_ptr<Meshes> meshes = nullptr;

  std::shared_ptr<CoarseGeometry> cGeo = nullptr;
public:
  explicit IProblem(std::string meshName = "") : meshesName(meshName) {
    if (meshesName.empty())
      Config::Get("Mesh", meshesName);

#ifndef BUILD_IA
    if (!meshesName.empty()) {
      // This can happen when problems are constructed and *later* initialized using
      // a MeshesCreator (spacetime).
      // TODO: Make those callers pass a MeshesCreator to a constructor instead.
      CreateMeshes(MeshesCreator(meshesName));
    }
#endif
    Config::Get("ProblemVerbose", verbose);
  }

  explicit IProblem(std::shared_ptr<Meshes> meshes) : meshes(std::move(meshes)) {
    Config::Get("ProblemVerbose", verbose);
  }

  explicit IProblem(IProblem &otherProblem) :
      meshes(otherProblem.meshes), verbose(otherProblem.verbose) {}

  virtual ~IProblem() = default;

  virtual bool HasExactSolution() const { return false; }

  void CreateMeshes(MeshesCreator meshesCreator) {
    std::string injectedMeshName = meshesCreator.GetMeshName();
    cGeo = CreateCoarseGeometryShared(injectedMeshName.empty() ? meshesName : injectedMeshName);
    meshes = meshesCreator.WithCoarseGeometry(cGeo).CreateShared();
  }

  const CoarseGeometry &Domain() const { return *cGeo; };

  const Meshes &GetMeshes() const {
    if (meshes == nullptr) Exit("Meshes not initialized") return *meshes;
  }

  bool IsMeshInitialized() const { return meshes != nullptr; }

  const Mesh &GetMesh(int level) const { return (*meshes)[level]; }

  const Mesh &GetMesh(LevelPair level) const { return (*meshes)[level]; }

  virtual std::string Name() const = 0;
};

#endif
