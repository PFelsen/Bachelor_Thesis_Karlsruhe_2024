We use the google benchmark library to run benchmarks:

https://github.com/google/benchmark


### Adding a benchmark
You can add a Benchmark by adding the cmake command

```add_mpp_bench(BenchFilename LIB_BENCH)```

where lib_bench is the needed library.

### Running a benchmark
First, call the benchmark you wish to perform by

```./MyBenchmark --benchmark_out=benchmark_filename --benchmark_out_format=json --benchmark_repetitions=10```

The repetitions assure you get a reasonable mean value. You can compare your benchmark to older versions with

```mpp/benchmark/tools/compare.py benchmarks old_benchmark_filename new_benchmark_filename```

Additionally, we should write a script evaluating the mean value differences...