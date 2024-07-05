#ifndef TESTENVIROMENT_HPP
#define TESTENVIROMENT_HPP

#include "Config.hpp"
#include "Parallel.hpp"
#include "Point.hpp"


#ifdef BUILD_IA

#include "IACInterval.hpp"


#endif

#include "benchmark/benchmark.h"
#include "gtest/gtest.h"

#include <complex>
#include <numeric>
#include <vector>


constexpr double tolerance = 1e-14;

void setSeed(int seed);

double RandomDouble(double min = -10, double max = 10);

double RandomNonZeroDouble(double min = -10, double max = 10);

int RandomInt(int min = -10, int max = 10);

int RandomNonZeroInt(int min = -10, int max = 10);

std::complex<double> RandomComplex(double min = -10, double max = 10);

std::complex<double> RandomNonZeroComplex(double min = -10, double max = 10);

using namespace ::testing;

typedef std::map<std::string, std::string> ConfigMap;

class MppTestBuilder {
public:
  int argc;

  char **argv;

  bool initPPM = false;

  std::string searchPath = "";
  std::string logPath = "";
  std::string geoPath = "";
  std::string plotPath = "";

  std::string configFile = "";

  bool initRandom = false;

  bool initDefaultConfig = true;

  bool parallelListeners = false;

  bool disableScreenLogging = true;

  bool disableFileLogging = true;


  ConfigMap defaultVerboseMap = {{"MeshVerbose", "1"},
                                 {"MeshesVerbose", "1"},
                                 {"ConfigVerbose", "0"},
                                 {"DistributionVerbose", "0"}};

  MppTestBuilder(int argc, char **argv) : argc(argc), argv(argv) {}

  MppTestBuilder WithPPM() {
    initPPM = true;
    return *this;
  }

  MppTestBuilder WithScreenLogging() {
    disableScreenLogging = false;
    return *this;
  }

  MppTestBuilder WithFileLogging() {
    disableFileLogging = false;
    return *this;
  }

  [[deprecated("Use WithConfPath instead")]]
  MppTestBuilder WithSearchPath(const std::string &_path_) {
    return WithConfPath(_path_);
  }

  MppTestBuilder WithConfPath(const std::string &_path_) {
    searchPath = _path_;
    return *this;
  }

  MppTestBuilder WithLogPath(const std::string &_path_) {
    logPath = _path_;
    return *this;
  }

  MppTestBuilder WithGeoPath(const std::string &_path_) {
    geoPath = _path_;
    return *this;
  }

  MppTestBuilder WithPlotPath(const std::string &_path_) {
    plotPath = _path_;
    return *this;
  }

  MppTestBuilder WithoutDefaultConfig() {
    initDefaultConfig = false;
    return *this;
  }

  MppTestBuilder WithConfigFile(const std::string &_confFile) {
    configFile = _confFile;
    initDefaultConfig = false;
    return *this;
  }

  MppTestBuilder WithConfigEntry(const std::string &key, const char *value) {
    initDefaultConfig = true;
    defaultVerboseMap[key] = std::string(value);
    return *this;
  }

  template<typename T>
  MppTestBuilder WithConfigEntry(const std::string &key, const T &value) {
    initDefaultConfig = true;
    defaultVerboseMap[key] = std::to_string(value);
    return *this;
  }

  MppTestBuilder WithParallelListeners() {
    parallelListeners = true;
    return *this;
  }

#ifdef BUILD_UQ
  MppTestBuilder WithRandomInitialized() {
    initRandom = true;
    return *this;
  }
#endif
};

class MppTest {
  bool withPPM;
public:
  MppTest(MppTestBuilder testBuilder) :
      MppTest(testBuilder.argc, testBuilder.argv, testBuilder.initPPM,
              testBuilder.disableScreenLogging, testBuilder.disableFileLogging,
              testBuilder.searchPath, testBuilder.logPath, testBuilder.geoPath,
              testBuilder.plotPath, testBuilder.configFile, testBuilder.initDefaultConfig,
              testBuilder.initRandom, testBuilder.parallelListeners,
              testBuilder.defaultVerboseMap) {}

  MppTest(int argc, char **argv, bool initPPM, bool disableScreenLogging, bool disableFileLogging,
          const std::string &searchPath, const std::string &logPath, const std::string &geoPath,
          const std::string &plotPath, const std::string &confFile, bool initDefaultConfig,
          bool initRandom, bool parallelListeners, ConfigMap defaultVerboseMap);

  int RUN_ALL_MPP_TESTS();

  int RUN_ALL_MPP_BENCHMARKS(int argc, char **argv);
};

inline void EXPECT_COMPLEX_EQ(const std::complex<double> &a, const std::complex<double> &b) {
  EXPECT_DOUBLE_EQ(std::real(a), std::real(b));
  EXPECT_DOUBLE_EQ(std::imag(a), std::imag(b));
}

void EXPECT_POINT_EQ(const Point &p1, const Point &p2);

void EXPECT_POINT_NE(const Point &p1, const Point &p2);

void EXPECT_POINT_NEAR(const Point &p1, const Point &p2, double tol);

template<typename T>
void EXPECT_VECTOR_EQ(const std::vector<T> &val1, const std::vector<T> &val2) {
  EXPECT_EQ(val1.size(), val2.size());
  for (int i = 0; i < val1.size(); ++i)
    EXPECT_EQ(val1[i], val2[i]);
}

template<typename T>
void EXPECT_VECTOR_NEAR(const std::vector<T> &val1, const std::vector<T> &val2, double tol) {
  EXPECT_EQ(val1.size(), val2.size());
  for (int i = 0; i < val1.size(); ++i)
    EXPECT_NEAR(val1[i], val2[i], tol) << " index: " << i;
}


#ifdef BUILD_IA

void EXPECT_IAINTERVAL_NEAR(const IAInterval &a, const IAInterval &b, double tol);

IAInterval RandomIAInterval(double min = -10, double max = 10, double widthMin = 1e-10,
                            double widthMax = 1e-1);

IAInterval RandomNonZeroIAInterval(double min = -10, double max = 10, double widthMin = 1e-10,
                                   double widthMax = 1e-1);

IACInterval RandomIACInterval(double min = -10, double max = 10, double widthMin = 1e-10,
                              double widthMax = 1e-1);

IACInterval RandomNonZeroIACInterval(double min = -10, double max = 10, double widthMin = 1e-10,
                                     double widthMax = 1e-1);

#endif

#endif // TESTENVIROMENT_HPP
