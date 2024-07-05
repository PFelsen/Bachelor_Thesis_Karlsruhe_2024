#include <random>
#include <MeshesCreator.hpp>
#include "benchmark/benchmark.h"

#include "TestEnvironment.hpp"

#include "Cell.hpp"

std::vector<Point> RandomCenters(int size) {
  std::uniform_real_distribution<double> unif(0, 5);
  std::default_random_engine re;

  std::vector<Point> centers(size);

  for (int i = 0; i < size; ++i) {
    centers.emplace_back(Point(unif(re), unif(re), unif(re)));
  }

  return centers;
}

static void BM_UnorderdMapUniqueCells(benchmark::State &state) {
  std::unordered_map<Point, std::unique_ptr<Cell>> meshCells{};

  int size = state.range(0);
  auto centers = RandomCenters(size);
  for (auto p : centers) {
    meshCells.try_emplace(p, std::unique_ptr<Cell>(
                                 CreateCell(QUADRILATERAL, 0,
                                            std::vector{p + Point(-1, -1, 0), p + Point(1, -1, 0),
                                                        p + Point(1, 1, 0), p + Point(-1, 1, 0)})));
  }

  // Random numbers to find random points

  std::uniform_int_distribution<int> unif(0, centers.size() - 1);
  std::default_random_engine re;

  // Perform setup here
  for (auto _ : state) {
    for (int i = 0; i < state.range(1); ++i) {
      auto p = centers[unif(re)];

      Cell &c = (*(meshCells.find(p))->second);
    }
  }
}

static void BM_UnorderdMapRawCells(benchmark::State &state) {
  std::unordered_map<Point, Cell *> meshCells{};

  int size = state.range(0);
  auto centers = RandomCenters(size);
  for (auto p : centers) {
    meshCells.try_emplace(p, CreateCell(QUADRILATERAL, 0,
                                        std::vector{p + Point(-1, -1, 0), p + Point(1, -1, 0),
                                                    p + Point(1, 1, 0), p + Point(-1, 1, 0)}));
  }

  // Random numbers to find random points

  std::uniform_int_distribution<int> unif(0, centers.size() - 1);
  std::default_random_engine re;

  // Perform setup here
  for (auto _ : state) {
    for (int i = 0; i < state.range(1); ++i) {
      auto p = centers[unif(re)];

      Cell &c = (*(meshCells.find(p))->second);
    }
  }
}

static void BM_UnorderdMapCellIterator(benchmark::State &state) {
  std::unordered_map<Point, Cell *> meshCells{};

  int size = state.range(0);
  auto centers = RandomCenters(size);
  for (auto p : centers) {
    meshCells.try_emplace(p, CreateCell(QUADRILATERAL, 0,
                                        std::vector{p + Point(-1, -1, 0), p + Point(1, -1, 0),
                                                    p + Point(1, 1, 0), p + Point(-1, 1, 0)}));
  }

  // Random numbers to find random points

  std::uniform_int_distribution<int> unif(0, centers.size() - 1);
  std::default_random_engine re;

  // Perform setup here
  for (auto _ : state) {
    for (int i = 0; i < state.range(1); ++i) {
      auto p = centers[unif(re)];

      cell c(meshCells.find(p));
    }
  }
}

static void BM_SortedVectorCells(benchmark::State &state) {
  std::vector<std::unique_ptr<Cell>> meshCells{};

  int size = state.range(0);
  auto centers = RandomCenters(size);
  for (auto p : centers) {
    meshCells.emplace_back(
        std::unique_ptr<Cell>(CreateCell(QUADRILATERAL, 0,
                                         std::vector{p + Point(-1, -1, 0), p + Point(1, -1, 0),
                                                     p + Point(1, 1, 0), p + Point(-1, 1, 0)})));
  }
  std::sort(meshCells.begin(), meshCells.end(),
            [](const std::unique_ptr<Cell> &a, const std::unique_ptr<Cell> &b) {
              return a->Center() < b->Center();
            });

  // Random numbers to find random points

  std::uniform_int_distribution<int> unif(0, centers.size() - 1);
  std::default_random_engine re;

  // Perform setup here
  for (auto _ : state) {
    for (int i = 0; i < state.range(1); ++i) {
      Point &p = centers[unif(re)];

      auto cellIt =
          std::find_if(meshCells.begin(), meshCells.end(),
                       [&p](const std::unique_ptr<Cell> &a) { return (a->Center()) == p; });
      Cell &c = *(meshCells[std::distance(meshCells.begin(), cellIt)]);
    }
  }
}

static void BM_UnsortedVectorCells(benchmark::State &state) {
  std::vector<std::unique_ptr<Cell>> meshCells{};

  int size = state.range(0);
  auto centers = RandomCenters(size);
  for (auto p : centers) {
    meshCells.emplace_back(
        std::unique_ptr<Cell>(CreateCell(QUADRILATERAL, 0,
                                         std::vector{p + Point(-1, -1, 0), p + Point(1, -1, 0),
                                                     p + Point(1, 1, 0), p + Point(-1, 1, 0)})));
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(meshCells.begin(), meshCells.end(), g);

  // Random numbers to find random points
  std::uniform_int_distribution<int> unif(0, centers.size() - 1);
  std::default_random_engine re;

  // Perform setup here
  for (auto _ : state) {
    for (int i = 0; i < state.range(1); ++i) {
      Point &p = centers[unif(re)];

      auto cellIt =
          std::find_if(meshCells.begin(), meshCells.end(),
                       [&p](const std::unique_ptr<Cell> &a) { return (a->Center()) == p; });
      Cell &c = *(meshCells[std::distance(meshCells.begin(), cellIt)]);
    }
  }
}

static void BM_MeshCells(benchmark::State &state) {
  auto meshes = MeshesCreator("Hexahedron").WithLevel(state.range(0)).CreateUnique();
  const Mesh &M = meshes->fine();


  std::vector<Point> centers(M.CellCount());
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    centers.emplace_back(c->first);
  }

  // Random numbers to find random points
  std::uniform_int_distribution<int> unif(0, centers.size() - 1);
  std::default_random_engine re;
}

BENCHMARK(BM_UnorderdMapUniqueCells)
    ->Ranges({{4 << 10, 16 << 10}, {16 << 10, 128 << 10}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_UnorderdMapRawCells)
    ->Ranges({{4 << 10, 16 << 10}, {16 << 10, 128 << 10}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_UnorderdMapCellIterator)
    ->Ranges({{4 << 10, 16 << 10}, {16 << 10, 128 << 10}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_SortedVectorCells)
    ->Ranges({{4 << 10, 16 << 10}, {2 << 10, 32 << 10}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_UnsortedVectorCells)
    ->Ranges({{4 << 10, 16 << 10}, {2 << 10, 32 << 10}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_MeshCells)
    ->ArgsProduct({{0, 1, 2, 3, 4, 5, 6}, {16 << 10, 32 << 10, 128 << 10}})
    ->Unit(benchmark::kMillisecond);

// Run the benchmark
int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_BENCHMARKS(argc, argv);
}