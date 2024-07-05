#include <random>
#include "VectorField.hpp"
#include "benchmark/benchmark.h"

#include "TestEnvironment.hpp"

static void BM_VectorFieldMultiplication(benchmark::State &state) {
  double lower_bound = 0;
  double upper_bound = 10000;
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  std::default_random_engine re;


  // Perform setup here
  for (auto _ : state) {
    state.PauseTiming();
    VectorField a(unif(re), unif(re), unif(re));
    VectorField b(unif(re), unif(re), unif(re));
    state.ResumeTiming();
    a *b;
  }
}

// Register the function as a benchmark
BENCHMARK(BM_VectorFieldMultiplication);

// Run the benchmark
int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_BENCHMARKS(argc, argv);
}
