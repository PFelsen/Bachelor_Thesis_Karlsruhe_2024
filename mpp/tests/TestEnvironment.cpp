#include "TestEnvironment.hpp"
#ifdef BUILD_UQ
#include "Random.hpp"
#endif
#include <VtuPlot.hpp>

void setSeed(int seed) { srand(seed); }

double RandomDouble(double min, double max) {
  return min + double(std::rand()) / RAND_MAX * (max - min);
}

double RandomNonZeroDouble(double min, double max) {
  double r = RandomDouble(min, max);
  while (r == 0.0) {
    r = RandomDouble(min, max);
  }
  return r;
}

int RandomInt(int min, int max) { return int(RandomDouble(int(min), int(max))); }

int RandomNonZeroInt(int min, int max) {
  int r = RandomInt(min, max);
  while (r == 0) {
    r = RandomInt(min, max);
  }
  return r;
}

std::complex<double> RandomComplex(double min, double max) {
  return std::complex<double>(RandomDouble(min, max), RandomDouble(min, max));
}

std::complex<double> RandomNonZeroComplex(double min, double max) {
  double re = RandomDouble(min, max);
  double im = RandomDouble(min, max);
  while (re * im == 0.0) {
    re = RandomDouble(min, max);
    im = RandomDouble(min, max);
  }
  return std::complex<double>(re, im);
}

MppTest::MppTest(int argc, char **argv, bool initPPM, bool disableScreenLogging,
                 bool disableFileLogging, const std::string &searchPath, const std::string &logPath,
                 const std::string &geoPath, const std::string &plotPath,
                 const std::string &confFile, bool initDefaultConfig, bool initRandom,
                 bool parallelListeners, ConfigMap defaultVerboseMap) {
  withPPM = initPPM;
  InitGoogleTest(&argc, argv);

  Config::SaveUsedConf(false);

  TestEventListeners &listeners = UnitTest::GetInstance()->listeners();
  if (initPPM) {
    ParallelProgrammingModel::Initialize(argc, argv);
    if (!PPM->Master(0) && !parallelListeners) {
      delete listeners.Release(listeners.default_result_printer());
    }
  }

  if (!searchPath.empty()) Config::SetConfPath(searchPath);
  if (!logPath.empty()) Config::SetLogPath(logPath);
  if (!geoPath.empty()) Config::SetGeoPath(geoPath);
  if (!plotPath.empty()) Config::SetPlotPath(plotPath);

  if (initDefaultConfig) {
    Config::Initialize(defaultVerboseMap);
  } else if (!confFile.empty()) {
    Config::Initialize(confFile);
    if (initPPM) {
      mout.setScreenEnabled(!disableScreenLogging);
      mout.setFileEnabled(!disableFileLogging);
    }
  }

  if (Config::IsInitialized()) {
    Config::Get("ParallelPlotting", globalParallelPlotting);
    Config::Get("Compression", globalCompression);
  }
#ifdef BUILD_UQ
  if (initRandom) { Random::Initialize(); }
#endif
}

int MppTest::RUN_ALL_MPP_TESTS() {
  int rc = RUN_ALL_TESTS();

  if (withPPM) ParallelProgrammingModel::Close();

  return rc;
}

int MppTest::RUN_ALL_MPP_BENCHMARKS(int argc, char **argv) {
  int status = -1;

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) status = 1;
  ::benchmark::RunSpecifiedBenchmarks();

  if (withPPM) ParallelProgrammingModel::Close();

  return status;
}

void EXPECT_POINT_EQ(const Point &p1, const Point &p2) {
  for (int i = 0; i < p1.SpaceDim() + p1.TimeDim(); ++i)
    EXPECT_EQ(p1[i], p2[i]);
}

void EXPECT_POINT_NE(const Point &p1, const Point &p2) {
  bool equal = true;
  for (int i = 0; i < p1.SpaceDim() + p1.TimeDim(); ++i) {
    equal = (p1[i] == p2[i]);
    if (!equal) break;
  }
  EXPECT_FALSE(equal);
}

void EXPECT_POINT_NEAR(const Point &p1, const Point &p2, double tol) {
  for (int i = 0; i < p1.SpaceDim() + p1.TimeDim(); ++i)
    EXPECT_NEAR(p1[i], p2[i], tol);
}

#ifdef BUILD_IA

void EXPECT_IAINTERVAL_NEAR(const IAInterval &a, const IAInterval &b, double tol) {
  EXPECT_NEAR(inf(a), inf(b), tol);
  EXPECT_NEAR(sup(a), sup(b), tol);
}

IAInterval RandomIAInterval(double min, double max, double widthMin, double widthMax) {
  double a = RandomDouble(min, max);
  return IAInterval(a, a + RandomDouble(widthMin, widthMax));
}

IAInterval RandomNonZeroIAInterval(double min, double max, double widthMin, double widthMax) {
  IAInterval r = RandomIAInterval(min, max, widthMin, widthMax);
  while (0.0 <= r) {
    r = RandomIAInterval(min, max, widthMin, widthMax);
  }
  return r;
}

IACInterval RandomIACInterval(double min, double max, double widthMin, double widthMax) {
  return IACInterval(RandomIAInterval(min, max, widthMin, widthMax),
                     RandomIAInterval(min, max, widthMin, widthMax));
}

IACInterval RandomNonZeroIACInterval(double min, double max, double widthMin, double widthMax) {
  IAInterval re = RandomIAInterval(min, max, widthMin, widthMax);
  IAInterval im = RandomIAInterval(min, max, widthMin, widthMax);
  while (0.0 <= re * im) {
    re = RandomIAInterval(min, max, widthMin, widthMax);
    im = RandomIAInterval(min, max, widthMin, widthMax);
  }
  return IACInterval(re, im);
}

#endif