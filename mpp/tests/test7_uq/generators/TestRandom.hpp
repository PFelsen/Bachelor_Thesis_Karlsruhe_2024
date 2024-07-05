#ifndef TESTUNIFORMDISTRIBUTION_HPP
#define TESTUNIFORMDISTRIBUTION_HPP

#include "MeshesCreator.hpp"
#include "Random.hpp"
#include "TestEnvironment.hpp"

class TestRandom : public Test {
protected:
  int numSamples = 1e6;

  int commSplit = 0;

  TestRandom(int commSplit) {
    if (commSplit > PPM->MaxCommSplit()) {
      this->commSplit = PPM->MaxCommSplit();
    } else {
      this->commSplit = commSplit;
    }
  }

  void TearDown() override { PPM->Barrier(0); }
};

class TestUniformDistribution : public TestRandom {
protected:
  double a;

  double b;

  template<typename T>
  void CheckResultsUniform(const T &results) const {
    EXPECT_NEAR((a + b) / 2.0, results.Mean(), sqrt(4.0 / numSamples));
    EXPECT_NEAR(pow(b - a, 2) / 12.0, results.Variance(), sqrt(10.0 / numSamples));
  };

  // Todo: implement statistical test

  TestUniformDistribution(int commSplit, double a, double b) : TestRandom(commSplit), a(a), b(b) {}

  void TearDown() override { PPM->Barrier(0); }
};

class TestUniformDistribution021 : public TestUniformDistribution {
public:
  explicit TestUniformDistribution021(int commSplit = 0) :
      TestUniformDistribution(commSplit, 0.0, 1.0) {}
};

class TestUniformDistribution021WithSplit : public TestUniformDistribution021 {
public:
  TestUniformDistribution021WithSplit() : TestUniformDistribution021(1) {}
};

class TestUniformDistribution021WithDoubleSplit : public TestUniformDistribution021 {
public:
  TestUniformDistribution021WithDoubleSplit() : TestUniformDistribution021(2) {}
};

class TestUniformDistribution021WithFullSplit : public TestUniformDistribution021 {
public:
  TestUniformDistribution021WithFullSplit() : TestUniformDistribution021(PPM->MaxCommSplit()) {}
};

class TestUniformDistributionNeg12Pos1 : public TestUniformDistribution {
public:
  explicit TestUniformDistributionNeg12Pos1(int commSplit = 0) :
      TestUniformDistribution(commSplit, -1.0, 1.0) {}
};

class TestUniformDistributionNeg12Pos1WithSplit : public TestUniformDistributionNeg12Pos1 {
public:
  TestUniformDistributionNeg12Pos1WithSplit() : TestUniformDistributionNeg12Pos1(1) {}
};

class TestUniformDistributionNeg12Pos1DoubleSplit : public TestUniformDistributionNeg12Pos1 {
public:
  TestUniformDistributionNeg12Pos1DoubleSplit() : TestUniformDistributionNeg12Pos1(2) {}
};

class TestUniformDistributionNeg12Pos1FullSplit : public TestUniformDistributionNeg12Pos1 {
public:
  TestUniformDistributionNeg12Pos1FullSplit() :
      TestUniformDistributionNeg12Pos1(PPM->MaxCommSplit()) {}
};

#define TEST_UNIFORM(TestClass)                                                                    \
                                                                                                   \
  TEST_F(TestClass, UniformDistributionReal) {                                                     \
    RVector results(numSamples);                                                                   \
    for (int i = 0; i < numSamples; i++)                                                           \
      results[i] = Random::Uniform(commSplit, a, b);                                               \
    CheckResultsUniform(results);                                                                  \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, UniformDistributionRVector) {                                                  \
    RVector results = Random::Uniform(commSplit, numSamples, a, b);                                \
    CheckResultsUniform(results);                                                                  \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, UniformDistributionRMatrix) {                                                  \
    RMatrix results = Random::Uniform(commSplit, sqrt(numSamples), sqrt(numSamples), a, b);        \
    CheckResultsUniform(results);                                                                  \
  }

class TestNormalDistribution : public TestRandom {
protected:
  double mean;

  double var;

  explicit TestNormalDistribution(int commSplit = 0, double mean = 0.0, double var = 1.0) :
      TestRandom(commSplit), mean(mean), var(var) {}

  // Todo: implement statistical test

  template<typename T>
  void CheckResultsNormal(const T &results) const {
      // TODO
  };
};

class TestNormalDistributionWithoutSplit : public TestNormalDistribution {
public:
  TestNormalDistributionWithoutSplit() : TestNormalDistribution(0) {}
};

class TestNormalDistributionWithSplit : public TestNormalDistribution {
public:
  TestNormalDistributionWithSplit() : TestNormalDistribution(1) {}
};

class TestNormalDistributionWithDoubleSplit : public TestNormalDistribution {
public:
  TestNormalDistributionWithDoubleSplit() : TestNormalDistribution(2) {}
};

class TestNormalDistributionWithFullSplit : public TestNormalDistribution {
public:
  TestNormalDistributionWithFullSplit() : TestNormalDistribution(PPM->MaxCommSplit()) {}
};

#define TEST_NORMAL(TestClass)                                                                     \
                                                                                                   \
  TEST_F(TestClass, NormalDistributionReal) {                                                      \
    RVector results(numSamples);                                                                   \
    for (int i = 0; i < numSamples; i++)                                                           \
      results[i] = Random::Normal(commSplit);                                                      \
    EXPECT_NEAR(0.0, results.Mean(), sqrt(10.0 / numSamples));                                     \
    EXPECT_NEAR(1.0, results.Variance(), sqrt(100.0 / numSamples));                                \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, ComplexNormalDistribution) {                                                   \
    CVector results(numSamples);                                                                   \
    for (int i = 0; i < numSamples; i++)                                                           \
      results[i] = Random::ComplexNormal(commSplit);                                               \
    Complex mean = results.Mean();                                                                 \
    EXPECT_NEAR(0.0, mean.real(), sqrt(10.0 / numSamples));                                        \
    EXPECT_NEAR(0.0, mean.imag(), sqrt(10.0 / numSamples));                                        \
    EXPECT_NEAR(2.0, results.Variance(), sqrt(100.0 / numSamples));                                \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, NormalDistributionRVector) {                                                   \
    RVector results = Random::Normal(commSplit, numSamples);                                       \
    EXPECT_NEAR(0.0, results.Mean(), sqrt(10.0 / numSamples));                                     \
    EXPECT_NEAR(1.0, results.Variance(), sqrt(100.0 / numSamples));                                \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, NormalDistributionCVector) {                                                   \
    CVector results = Random::ComplexNormal(commSplit, numSamples);                                \
    Complex mean = results.Mean();                                                                 \
    EXPECT_NEAR(0.0, mean.real(), sqrt(100.0 / numSamples));                                       \
    EXPECT_NEAR(0.0, mean.imag(), sqrt(100.0 / numSamples));                                       \
    EXPECT_NEAR(2.0, results.Variance(), sqrt(100.0 / numSamples));                                \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, NormalDistributionRMatrix) {                                                   \
    RMatrix results = Random::Normal(commSplit, int(sqrt(numSamples)), int(sqrt(numSamples)));     \
    EXPECT_NEAR(0.0, results.Mean(), sqrt(100.0 / numSamples));                                    \
    EXPECT_NEAR(1.0, results.Variance(), sqrt(100.0 / numSamples));                                \
  }                                                                                                \
                                                                                                   \
  TEST_F(TestClass, NormalDistributionCMatrix) {                                                   \
    CMatrix results =                                                                              \
        Random::ComplexNormal(commSplit, int(sqrt(numSamples)), int(sqrt(numSamples)));            \
    Complex mean = results.Mean();                                                                 \
    EXPECT_NEAR(0.0, mean.real(), sqrt(100.0 / numSamples));                                       \
    EXPECT_NEAR(0.0, mean.imag(), sqrt(100.0 / numSamples));                                       \
    EXPECT_NEAR(2.0, results.Variance(), sqrt(100.0 / numSamples));                                \
  }

// TEST_F(TestClass, NormalDistributionRTensor) {
//   RTensor results =
//       Random::Normal(commSplit, int(cbrt(numSamples)),
//                      int(cbrt(numSamples)), int(cbrt(numSamples)));
//   EXPECT_NEAR(0.0, results.Mean(), sqrt(100.0 / numSamples));
//   EXPECT_NEAR(1.0, results.Variance(), sqrt(100.0 / numSamples));
// }
//
// TEST_F(TestClass, NormalDistributionCTensor) {
//   CTensor results =
//       Random::ComplexNormal(commSplit, int(cbrt(numSamples)),
//                             int(cbrt(numSamples)), int(cbrt(numSamples)));
//   Complex mean = results.Mean();
//   EXPECT_NEAR(0.0, mean.real(), sqrt(100.0 / numSamples));
//   EXPECT_NEAR(0.0, mean.imag(), sqrt(100.0 / numSamples));
//   EXPECT_NEAR(2.0, results.Variance(), sqrt(100.0 / numSamples));


#endif // TESTUNIFORMDISTRIBUTION_HPP
