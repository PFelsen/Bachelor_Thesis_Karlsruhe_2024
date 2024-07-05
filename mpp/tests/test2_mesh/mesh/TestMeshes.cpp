#include <cstddef>
#include <memory>
#include <string>
#include <utility>

#include <CoarseGeometry.hpp>
#include <Distribution.hpp>
#include <LevelPair.hpp>
#include <Mesh.hpp>
#include <MeshSettings.hpp>
#include <Meshes.hpp>
#include <MeshesCreator.hpp>
#include <Parallel.hpp>
#include <gtest/gtest.h>

#include "TestEnvironment.hpp"
#include "TestMeshesPrintInfo.hpp"

class TestMeshes : public Meshes {
public:
  explicit TestMeshes(MeshSettings &&settings) : Meshes(std::move(settings)) {}

  MeshesContainer::node_type extract(const LevelPair &lp) { return meshes.extract(lp); }

  void insert(MeshesContainer::node_type &&node) { meshes.insert(std::move(node)); }
};

TEST(TestMeshes, TestSpaceTimeMeshes) {
#ifndef USE_SPACETIME
  GTEST_SKIP() << "SpaceTime is not enabled";
#endif

  static constexpr int scales = 2;
  static constexpr int coarseLevel = 3;
  static constexpr int distributeLevel = 2;
  static constexpr int fineLevel = 4;

  static constexpr LevelPair coarseLp = LevelPair{coarseLevel, -1};
  static constexpr LevelPair distributedLp = LevelPair{distributeLevel, distributeLevel};
  static constexpr LevelPair fineLp = LevelPair{fineLevel, fineLevel};
  static constexpr LevelPair maxDefaultLp = LevelPair{fineLevel, fineLevel + scales};

  MeshSettings settings;
  settings.coarseLevel = coarseLevel;
  settings.distributeLevel = distributeLevel;
  settings.fineLevel = fineLevel;
  settings.timeRefinement = scales;
  settings.coarseGeometry =
      CreateCoarseGeometryShared("SpaceTimeSquare", VertexDataList(settings.vertexData),
                                 CellDataList(settings.cellData));

  TestMeshes meshes(std::move(settings));

  // Test if all meshes are present
  EXPECT_EQ(meshes.Size(), (fineLevel - distributeLevel + 1) /* Space only meshes */
                               + (fineLevel - distributeLevel + 1)
                                     * (fineLevel - distributeLevel + 1 + scales) /* ST meshes */);
  EXPECT_TRUE(meshes.Contains(coarseLp));
  EXPECT_TRUE(meshes.Contains(LevelPair{fineLevel, -1}));
  EXPECT_TRUE(meshes.Contains(distributedLp));
  EXPECT_TRUE(meshes.Contains(LevelPair{fineLevel, distributeLevel}));
  EXPECT_TRUE(meshes.Contains(LevelPair{distributeLevel, fineLevel}));
  EXPECT_TRUE(meshes.Contains(fineLp));
  EXPECT_TRUE(meshes.Contains(LevelPair{distributeLevel, fineLevel + scales}));
  EXPECT_TRUE(meshes.Contains(maxDefaultLp));

  // Test lazy init for > fineLevel & erase
  const std::size_t defaultSize = meshes.Size();
  static constexpr LevelPair beyondMaxDefaultLp = LevelPair{fineLevel + 1, fineLevel + scales + 1};
  {
    const Mesh &newMesh = meshes[beyondMaxDefaultLp];
    EXPECT_EQ(newMesh.Level(), beyondMaxDefaultLp);
  }
  /* Added: {fineLevel + 1, fineLevel + scales},
            {fineLevel, fineLevel + scales + 1},
            {fineLevel + 1, fineLevel + scales + 1} */
  EXPECT_EQ(meshes.Size(), defaultSize + 3);
  EXPECT_EQ(meshes.Erase(beyondMaxDefaultLp), 1);
  EXPECT_EQ(meshes.Size(), defaultSize + 2);

  // Test if meshes are distributed
  EXPECT_EQ((meshes[LevelPair{distributeLevel, -1}].ProcSetsCount()), 0);
  EXPECT_GT((meshes[distributedLp].ProcSetsCount()), 0);
  EXPECT_GT((meshes[fineLp].ProcSetsCount()), 0);
}

TEST(TestMeshes, TestSpaceMeshes) {
  static constexpr int coarseLevel = 3;
  static constexpr int distributeLevel = 2;
  static constexpr int fineLevel = 4;

  static constexpr LevelPair coarseLp = LevelPair{coarseLevel, -1};
  static constexpr LevelPair distributeLp = LevelPair{distributeLevel, -1};
  static constexpr LevelPair fineLp = LevelPair{fineLevel, -1};

  MeshSettings settings;
  settings.coarseLevel = coarseLevel;
  settings.distributeLevel = distributeLevel;
  settings.fineLevel = fineLevel;
  settings.coarseGeometry =
      CreateCoarseGeometryShared("Interval", VertexDataList(settings.vertexData),
                                 CellDataList(settings.cellData));

  TestMeshes meshes(std::move(settings));

  // Test if all meshes are present
  const std::size_t defaultSize = meshes.Size();
  EXPECT_EQ(defaultSize, fineLevel - distributeLevel + 1);
  EXPECT_TRUE(meshes.Contains(coarseLp));
  EXPECT_TRUE(meshes.Contains(distributeLp));
  EXPECT_TRUE(meshes.Contains(fineLp));

  // Test lazy init for > fineLevel & erase
  static constexpr LevelPair beyondFineLp = LevelPair{fineLevel + 1, -1, 0, 0};
  {
    const Mesh &newMesh = meshes[beyondFineLp];
    EXPECT_EQ(newMesh.Level(), beyondFineLp);
  }
  EXPECT_EQ(meshes.Size(), defaultSize + 1);
  EXPECT_EQ(meshes.Erase(beyondFineLp), 1);
  EXPECT_EQ(meshes.Size(), defaultSize);

  // Test init for new commSplits
  PPM->SplitCommunicators(1);
  PPM->Barrier(1);
  static constexpr LevelPair newCommSplitLevel = LevelPair{coarseLevel, -1, 0, 1};
  const Mesh &newCommSplitMesh = meshes[newCommSplitLevel];
  EXPECT_EQ(newCommSplitMesh.Level(), newCommSplitLevel);
  EXPECT_EQ(meshes.Size(), 2 * defaultSize);
  EXPECT_NE(&newCommSplitMesh, &(meshes[LevelPair{coarseLevel, -1, 0, 0}]));
  PPM->ClearCommunicators();
  PPM->Barrier(0);

  // Test ignore adaptivityLevel
  const Mesh *const meshPtrWithoutAdaptivity = &meshes[LevelPair{fineLevel, -1, 0, 0}];
  const Mesh *const meshPtrWithAdaptivity = &meshes[LevelPair{fineLevel, -1, 1, 0}];
  EXPECT_EQ(meshPtrWithoutAdaptivity, meshPtrWithAdaptivity);

  // Test if meshes are distributed
  EXPECT_GT((meshes[distributeLp].ProcSetsCount()), 0);
  EXPECT_GT((meshes[fineLp].ProcSetsCount()), 0);

  // Test distribution of lazy created meshes
  auto coarseNode = meshes.extract(coarseLp);
  auto distributedNode = meshes.extract(distributeLp);
  meshes.Clear();
  meshes.insert(std::move(coarseNode));
  const auto &newDistributedMesh = meshes[LevelPair{distributeLevel, distributeLevel}];
  EXPECT_EQ(distributedNode.mapped().Level(), newDistributedMesh.Level());
  EXPECT_EQ(distributedNode.mapped().CellCount(), newDistributedMesh.CellCount());
  EXPECT_EQ(distributedNode.mapped().VertexCount(), newDistributedMesh.VertexCount());
  EXPECT_EQ(distributedNode.mapped().FaceCount(), newDistributedMesh.FaceCount());
  EXPECT_EQ(distributedNode.mapped().EdgeCount(), newDistributedMesh.EdgeCount());
  EXPECT_EQ(distributedNode.mapped().BoundaryFaceCount(), newDistributedMesh.BoundaryFaceCount());
  EXPECT_EQ(distributedNode.mapped().ProcSetsCount(), newDistributedMesh.ProcSetsCount());
}

TEST(TestMeshes, TestSingleMeshMeshes) {
  static constexpr int level = 0;
  static constexpr LevelPair lp = LevelPair{level, -1};

  MeshSettings settings;
  settings.coarseLevel = level;
  settings.distributeLevel = level;
  settings.fineLevel = level;
  settings.coarseGeometry =
      CreateCoarseGeometryShared("Interval", VertexDataList(settings.vertexData),
                                 CellDataList(settings.cellData));

  TestMeshes meshes(std::move(settings));
  EXPECT_EQ(meshes.Size(), level + 1);
  EXPECT_TRUE(meshes.Contains(lp));
  EXPECT_GT((meshes[lp].ProcSetsCount()), 0);
}

TEST_P(TestMeshesPrintInfo, TestMeshesPrintInfo) { meshes->PrintInfo(); }

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("DistributionVerbose", 1)
                     .WithConfigEntry("MeshesVerbose", 4)
                     .WithConfigEntry("MeshVerbose", 2)
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
