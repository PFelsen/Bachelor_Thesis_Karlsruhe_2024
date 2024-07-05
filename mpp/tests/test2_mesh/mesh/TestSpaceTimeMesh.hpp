#ifndef SPACETIME_TESTSPACETIMEMESH_HPP
#define SPACETIME_TESTSPACETIMEMESH_HPP

#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"

using NumberOfProc = int;
using OnProc = int;

struct MeshDistributionOnProc {
  vector<Point> cell_points;
  vector<Point> overlap_points;
};

struct MeshDistribution {
  std::map<OnProc, MeshDistributionOnProc> per_proc;
};

MeshDistributionOnProc createMeshDistributionOnProc0for1() {
  vector<Point> cell_points{Point(0.25, 0.25, 0.125), Point(0.75, 0.25, 0.125),
                            Point(0.25, 0.75, 0.125), Point(0.75, 0.75, 0.125),
                            Point(0.25, 0.25, 0.375), Point(0.75, 0.25, 0.375),
                            Point(0.25, 0.75, 0.375), Point(0.75, 0.75, 0.375),
                            Point(0.25, 0.25, 0.625), Point(0.75, 0.25, 0.625),
                            Point(0.25, 0.75, 0.625), Point(0.75, 0.75, 0.625),
                            Point(0.25, 0.25, 0.875), Point(0.75, 0.25, 0.875),
                            Point(0.25, 0.75, 0.875), Point(0.75, 0.75, 0.875)};

  vector<Point> overlap_points{};

  return MeshDistributionOnProc{cell_points, overlap_points};
}

MeshDistributionOnProc createMeshDistributionOnProc0for2() {
  vector<Point> cell_points{Point(0.25, 0.25, 0.125), Point(0.75, 0.25, 0.125),
                            Point(0.25, 0.75, 0.125), Point(0.75, 0.75, 0.125),
                            Point(0.25, 0.25, 0.375), Point(0.75, 0.25, 0.375),
                            Point(0.25, 0.75, 0.375), Point(0.75, 0.75, 0.375)};

  vector<Point> overlap_points{Point(0.25, 0.25, 0.625), Point(0.75, 0.25, 0.625),
                               Point(0.25, 0.75, 0.625), Point(0.75, 0.75, 0.625)};

  return MeshDistributionOnProc{cell_points, overlap_points};
}

MeshDistributionOnProc createMeshDistributionOnProc1for2() {
  vector<Point> cell_points{Point(0.25, 0.25, 0.625), Point(0.75, 0.25, 0.625),
                            Point(0.25, 0.75, 0.625), Point(0.75, 0.75, 0.625),
                            Point(0.25, 0.25, 0.875), Point(0.75, 0.25, 0.875),
                            Point(0.25, 0.75, 0.875), Point(0.75, 0.75, 0.875)};
  vector<Point> overlap_points{Point(0.25, 0.25, 0.375), Point(0.75, 0.25, 0.375),
                               Point(0.25, 0.75, 0.375), Point(0.75, 0.75, 0.375)};
  return MeshDistributionOnProc{cell_points, overlap_points};
}

MeshDistribution createMeshDistribution1() {
  std::map<OnProc, MeshDistributionOnProc> per_proc;
  per_proc.insert({0, createMeshDistributionOnProc0for1()});
  return {per_proc};
}

MeshDistribution createMeshDistribution2() {
  std::map<OnProc, MeshDistributionOnProc> per_proc;
  per_proc.insert({0, createMeshDistributionOnProc0for2()});
  per_proc.insert({1, createMeshDistributionOnProc1for2()});
  return {per_proc};
}

std::map<NumberOfProc, MeshDistribution> createDisc() {
  std::map<NumberOfProc, MeshDistribution> dist;
  dist.insert({1, createMeshDistribution1()});
  dist.insert({2, createMeshDistribution2()});
  return dist;
}

struct SpaceTimeMeshTestParam {
  std::string name;
  std::map<NumberOfProc, MeshDistribution> distribution_data;
};

class TestSpaceTimeMesh : public TestWithParam<SpaceTimeMeshTestParam> {
public:
  std::string name;
  std::unique_ptr<Meshes> meshes;
  std::map<NumberOfProc, MeshDistribution> distribution_data;
  MeshDistributionOnProc dist_on_proc;

  TestSpaceTimeMesh() : name(GetParam().name), distribution_data(GetParam().distribution_data) {
    meshes = MeshesCreator(name).WithPLevel(1).WithLevel(1).CreateUnique();
    dist_on_proc = distribution_data[PPM->Size()].per_proc[PPM->Proc()];
  }

  void checkCellPoints() {
    const Mesh &mesh = (*meshes)[{1, 1}];
    EXPECT_EQ(dist_on_proc.cell_points.size(), mesh.CellCount());
    for (Point &p : dist_on_proc.cell_points) {
      EXPECT_TRUE(mesh.find_cell(p) != mesh.cells_end());
    }
  }

  void checkOverlapPoints() {
    const Mesh &mesh = (*meshes)[{1, 1}];
    EXPECT_EQ(dist_on_proc.overlap_points.size(), mesh.OverlapCount());
    for (Point &p : dist_on_proc.overlap_points) {
      EXPECT_TRUE(mesh.find_overlap_cell(p) != mesh.overlap_end());
    }
  }
};

INSTANTIATE_TEST_SUITE_P(TestSpaceTimeMesh, TestSpaceTimeMesh,
                         Values(SpaceTimeMeshTestParam{"SpaceTimeSquare", createDisc()}));


#endif // SPACETIME_TESTSPACETIMEMESH_HPP
