#ifndef TESTSTOCHASTICCOLLOCATION_HPP
#define TESTSTOCHASTICCOLLOCATION_HPP

#include "SingleLevelEstimator.hpp"
#include "TestEnvironment.hpp"


constexpr double SC_TEST_TOLERANCE = 1e-6;

struct TestParams {
  std::string problemName;

  std::string quantity;

  std::string model;

  double refValue;

  int level = 3;

  bool onlyFine = true;
};

Logging &operator<<(Logging &s, const TestParams &testParams) {
  return s << "Problem Name: " << testParams.problemName << endl
           << "Quantity: " << testParams.quantity << endl
           << "Model: " << testParams.model << endl;
}

class TestStochasticCollocation : public TestWithParam<TestParams> {
protected:
  int stochLevel;

  double epsilon;

  PDESolverConfig conf;

  std::unique_ptr<Estimator> scSeriell;

  std::unique_ptr<Estimator> scParallel;

  TestStochasticCollocation(double epsilon, int stochLevel) :
      stochLevel(stochLevel), epsilon(epsilon),

      conf(PDESolverConfig(GetParam().level, 1, GetParam().quantity, "size", GetParam().model,
                           GetParam().problemName)),


      scSeriell(EstimatorCreator()
                    .WithEstimator(STOCHASTIC_COLLOCATION)
                    .WithOnlyFine(GetParam().onlyFine)
                    .WithInitLevel(GetParam().level)
                    .WithPDESolverConfig(conf)
                    .WithStochLevel(stochLevel)
                    .WithEpsilon(epsilon)
                    .WithParallel(false)
                    .CreateUnique()),

      scParallel(EstimatorCreator()
                     .WithEstimator(STOCHASTIC_COLLOCATION)
                     .WithPDESolverConfig(conf)
                     .WithOnlyFine(GetParam().onlyFine)
                     .WithInitLevel(GetParam().level)
                     .WithStochLevel(stochLevel)
                     .WithEpsilon(epsilon)
                     .WithParallel(true)
                     .CreateUnique()) {}

  void TearDown() override { PPM->Barrier(0); }
};

class TestStochasticCollocationWithoutEpsilon : public TestStochasticCollocation {
public:
  TestStochasticCollocationWithoutEpsilon() : TestStochasticCollocation(0.0, 6) {}
};

class TestStochasticCollocationFailingTests : public TestStochasticCollocation {
public:
  TestStochasticCollocationFailingTests() : TestStochasticCollocation(0.0, 6) {}
};

#endif // TESTSTOCHASTICCOLLOCATION_HPP
