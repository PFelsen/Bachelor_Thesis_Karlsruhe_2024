#ifndef TESTSPARSEGRIDGENERATOR_HPP
#define TESTSPARSEGRIDGENERATOR_HPP

#include "MeshesCreator.hpp"
#include "SparseGridGenerator.hpp"
#include "TestEnvironment.hpp"


constexpr double TEST_TOLERANCE = 1e-10;

class TestSparseGridGenerator : public TestWithParam<int> {
protected:
  std::shared_ptr<Meshes> meshes;

  SparseGridGenerator generator;

  explicit TestSparseGridGenerator(const GridDomain &domain, bool local = false) :
      meshes(MeshesCreator("Interval").CreateShared()),
      generator(SparseGridGenerator(domain, GetParam(), 0, local)) {
    generator.DrawSample(SampleID(GetParam(), 0, false));
  }

  void TearDown() override {}
};

class TestClenshawCurtis1D : public TestSparseGridGenerator {
public:
  TestClenshawCurtis1D() : TestSparseGridGenerator(GridDomain(1)) {}
};

class TestClenshawCurtis2D : public TestSparseGridGenerator {
public:
  TestClenshawCurtis2D() : TestSparseGridGenerator(GridDomain(2)) {}
};

class TestClenshawCurtis5D : public TestSparseGridGenerator {
public:
  TestClenshawCurtis5D() : TestSparseGridGenerator(GridDomain(5)) {}
};

class TestClenshawCurtis1221D : public TestSparseGridGenerator {
public:
  TestClenshawCurtis1221D() : TestSparseGridGenerator(GridDomain({{1}, {2}})) {}
};

class TestClenshawCurtis0212D : public TestSparseGridGenerator {
public:
  TestClenshawCurtis0212D() : TestSparseGridGenerator(GridDomain({{0, 0}, {1, 1}})) {}
};

class TestGaussHermite1D : public TestSparseGridGenerator {
public:
  TestGaussHermite1D() : TestSparseGridGenerator(GridDomain(1, TasGrid::rule_gausshermite)) {}
};

class TestLocalPolynomial1D : public TestSparseGridGenerator {
public:
  TestLocalPolynomial1D() : TestSparseGridGenerator(GridDomain(1), true) {}
};

class TestLocalPolynomial2D : public TestSparseGridGenerator {
public:
  TestLocalPolynomial2D() : TestSparseGridGenerator(GridDomain(2), true) {}
};

#define SPARSEGRIDGENERATOR(TestClass)                                                             \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestSparseGridGenerator, TestClass, Values(0, 1, 2, 3, 4, 5, 6, 7, 8)); \
                                                                                                   \
  TEST_P(TestClass, TestNumOfGridPoints) {                                                         \
    EXPECT_EQ(generator.GetNumPoints(), generator.GetWeights().size());                            \
    EXPECT_EQ(generator.GetNumPoints(), generator.GetSamples().size());                            \
  }                                                                                                \
                                                                                                   \
  TEST_P(TestClass, TestSumOfWeights) {                                                            \
    double sum = 0.0;                                                                              \
    for (auto &weight : generator.GetWeights())                                                    \
      sum += weight;                                                                               \
    EXPECT_NEAR(sum, generator.SumOfWeights(), TEST_TOLERANCE);                                    \
  }                                                                                                \
                                                                                                   \
  TEST_P(TestClass, TestLinear) {                                                                  \
    RVector sum(generator.Domain().dim);                                                           \
    std::vector<double> weights = generator.GetWeights();                                          \
    std::vector<RVector> samples = generator.GetSamples();                                         \
    double n = generator.GetNumPoints();                                                           \
    for (int i = 0; i < n; i++) {                                                                  \
      sum += weights[i] * samples[i];                                                              \
    }                                                                                              \
    for (int d = 0; d < generator.Domain().dim; d++) {                                             \
      double value =                                                                               \
          0.5 * (pow(generator.Domain().right()[d], 2) - pow(generator.Domain().left()[d], 2));    \
      EXPECT_NEAR(sum[d], value, TEST_TOLERANCE);                                                  \
    }                                                                                              \
  }                                                                                                \
                                                                                                   \
  TEST_P(TestClass, TestGridPointsInDomain) {                                                      \
    GridDomain domain = generator.Domain();                                                        \
    for (auto &sample : generator.GetSamples()) {                                                  \
      for (int d = 0; d < domain.dim; d++) {                                                       \
        EXPECT_LE(domain.left()[d], sample[d]);                                                    \
        EXPECT_LE(sample[d], domain.right()[d]);                                                   \
      }                                                                                            \
    }                                                                                              \
  }


#endif // TESTSPARSEGRIDGENERATOR_HPP
