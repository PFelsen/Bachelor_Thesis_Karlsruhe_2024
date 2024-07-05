#include "TestEnvironment.hpp"

#include <unistd.h>
#include "LagrangeDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "MultiSampleFEM.hpp"
#include "PDESolver.hpp"
#include "SingleLevelEstimator.hpp"
#include "WelfordAggregate.hpp"

/*
 * Todo:
 *  - Design meaningful test for Cost Estimation
 *    -> EXPECT_NEAR(estimator.GetAggregate().mean.C : estimator.C, referenceMean, 0.01);
 *    -> Use sleep(PPM->Proc(0)); in PDESolver.run
 *  - Design meaningful test for Index also within while loop.
 */


class TestSLEstimator : public TestWithParam<int> {
protected:
  SingleLevelEstimator estimator;

  explicit TestSLEstimator(SLEstimatorConfig slEstmConfig) :
      estimator(SingleLevelEstimator(slEstmConfig.WithInitSamples(GetParam()))) {}

  void TestSVarExample1() {
    pout << "Before 1. Run: " << estimator.GetAggregate();
    estimator.Method();
    pout << "After 1. Run: " << estimator.GetAggregate();

    double referenceMean = 0.0;
    double referenceSVar = 0.0;
    for (int i = 0; i < GetParam(); i++) {
      if (i % 2) referenceSVar += pow(1.0 - referenceMean, 2);
      else referenceSVar += pow(-1.0 - referenceMean, 2);
    }
    referenceSVar = referenceSVar / (GetParam() - 1);

    EXPECT_EQ(estimator.GetAggregate().Q.GetSVar(), referenceSVar);
    EXPECT_EQ(estimator.GetAggregate().Q.GetMean(), referenceMean);
  }

  int ReferenceIndex(int start, int samples) {
    if (samples >= PPM->Size(0)) return start + (estimator.SamplesToComputeOnComm()) * PPM->Proc(0);
    else return start + PPM->Color(estimator.CommSplit());
  }

  void TestIndices(int start, int samples) {
    int referenceIndex = ReferenceIndex(start, samples);
    pout << "Before 1. Run: " << DOUT(estimator.Index()) << endl;
    EXPECT_EQ(estimator.Index(), referenceIndex);
    estimator.Method();
    pout << "After 1. Run: " << DOUT(estimator.Index()) << endl;
  }

  void TearDown() override { PPM->Barrier(0); }
};

class TestSLEstmProcEvaluation : public TestSLEstimator {
public:
  TestSLEstmProcEvaluation() :
      TestSLEstimator(SLEstimatorConfig().WithInitLevel(2).WithMultiSampleFEMConfig(
          MultiSampleFEMConfig("EmptyLagrangeProtocol", "ProcEvaluation")
              .WithPDESolverConfig(PDESolverConfig().WithDegree(1)))) {}
};

TEST_P(TestSLEstmProcEvaluation, TestInstanciation) {
  switch (GetParam()) {
  case 4:
    switch (PPM->Size(0)) {
    case 8:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 1);
      EXPECT_EQ(estimator.CommSplit(), 2);
      break;
    case 4:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 1);
      EXPECT_EQ(estimator.CommSplit(), 2);
      break;
    case 2:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 2);
      EXPECT_EQ(estimator.CommSplit(), 1);
      break;
    case 1:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 4);
      EXPECT_EQ(estimator.CommSplit(), 0);
      break;
    default:
      break;
    }
    break;
  case 8:
    switch (PPM->Size(0)) {
    case 8:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 1);
      EXPECT_EQ(estimator.CommSplit(), 3);
      break;
    case 4:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 2);
      EXPECT_EQ(estimator.CommSplit(), 2);
      break;
    case 2:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 4);
      EXPECT_EQ(estimator.CommSplit(), 1);
      break;
    case 1:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 8);
      EXPECT_EQ(estimator.CommSplit(), 0);
      break;
    default:
      break;
    }
    break;
  case 16:
    switch (PPM->Size(0)) {
    case 8:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 2);
      EXPECT_EQ(estimator.CommSplit(), 3);
      break;
    case 4:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 4);
      EXPECT_EQ(estimator.CommSplit(), 2);
      break;
    case 2:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 8);
      EXPECT_EQ(estimator.CommSplit(), 1);
      break;
    case 1:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 16);
      EXPECT_EQ(estimator.CommSplit(), 0);
      break;
    default:
      Warning("No test case for this amount of processes") break;
    }
    break;
  case 1000:
    switch (PPM->Size(0)) {
    case 8:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 125);
      EXPECT_EQ(estimator.CommSplit(), 3);
      break;
    case 4:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 250);
      EXPECT_EQ(estimator.CommSplit(), 2);
      break;
    case 2:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 500);
      EXPECT_EQ(estimator.CommSplit(), 1);
      break;
    case 1:
      pout << DOUT(estimator.GetAggregate().ctr) << DOUT(estimator.CommSplit()) << endl;
      EXPECT_EQ(estimator.SamplesToComputeOnComm(), 1000);
      EXPECT_EQ(estimator.CommSplit(), 0);
      break;
    default:
      break;
    }
    break;
  }
}

TEST_P(TestSLEstmProcEvaluation, TestMean) {
  if (PPM->Size(0) == 8 && GetParam() == 4) { return; } // The test is not meaningful for this case?

  estimator.Method();
  pout << "After 1. Run: " << estimator.GetAggregate().Q.GetTotal();

  double sumGauss = (pow((PPM->Size(0) - 1), 2) + (PPM->Size(0) - 1)) / 2.0;
  double referenceMean = sumGauss / PPM->Size(0);

  EXPECT_EQ(estimator.GetAggregate().Q.GetMean(), referenceMean);

  int additionalSamples = PPM->Size(0);
  estimator.UpdateSampleAmount(additionalSamples);

  pout << "Before 2. Run: " << estimator.GetAggregate().Q.GetMean();
  estimator.Method();


  pout << "After 2. Run: " << estimator.GetAggregate().Q.GetMean();

  EXPECT_EQ(estimator.GetAggregate().Q.GetMean(), referenceMean);
}

class TestSLEstmZeroEvaluation : public TestSLEstimator {
public:
  TestSLEstmZeroEvaluation() :
      TestSLEstimator(SLEstimatorConfig().WithInitLevel(2).WithMultiSampleFEMConfig(
          MultiSampleFEMConfig("EmptyLagrangeProtocol", "ZeroEvaluation")
              .WithPDESolverConfig(PDESolverConfig().WithDegree(1)))) {}
};

TEST_P(TestSLEstmZeroEvaluation, TestIndex) {
  int totalSamples = GetParam();
  pout << "Start with " << GetParam() << " samples" << endl;
  TestIndices(0, GetParam());
  EXPECT_EQ(totalSamples, estimator.GetAggregate().ctr.M);

  int additionalSamples = 4;
  totalSamples += additionalSamples;
  estimator.UpdateSampleAmount(additionalSamples);
  pout << "Add " << additionalSamples << " more samples" << endl;
  TestIndices(estimator.GetAggregate().ctr.M, additionalSamples);
  EXPECT_EQ(totalSamples, estimator.GetAggregate().ctr.M);

  additionalSamples = 8;
  totalSamples += additionalSamples;
  estimator.UpdateSampleAmount(additionalSamples);
  pout << "Add " << additionalSamples << " more samples" << endl;
  TestIndices(estimator.GetAggregate().ctr.M, additionalSamples);
  EXPECT_EQ(totalSamples, estimator.GetAggregate().ctr.M);

  additionalSamples = 16;
  totalSamples += additionalSamples;
  estimator.UpdateSampleAmount(additionalSamples);
  pout << "Add " << additionalSamples << " more samples" << endl;
  TestIndices(estimator.GetAggregate().ctr.M, additionalSamples);
  EXPECT_EQ(totalSamples, estimator.GetAggregate().ctr.M);
  mout << endl;
}

TEST_P(TestSLEstmZeroEvaluation, TestSampleCounter) {
  pout << "Before: " << estimator.GetAggregate().ctr;
  estimator.Method();
  pout << "After: " << estimator.GetAggregate().ctr;
  EXPECT_EQ(estimator.GetAggregate().ctr.M, GetParam());
  mout << endl;
}

class TestSLEstmModuloIndexEvaluation : public TestSLEstimator {
public:
  TestSLEstmModuloIndexEvaluation() :
      TestSLEstimator(SLEstimatorConfig().WithInitLevel(2).WithMultiSampleFEMConfig(
          MultiSampleFEMConfig("EmptyLagrangeProtocol", "ModuloIndexEvaluation")
              .WithPDESolverConfig(PDESolverConfig().WithDegree(1)))) {}
};

TEST_P(TestSLEstmModuloIndexEvaluation, TestSVar) {
  switch (GetParam()) {
  case 4:
    switch (PPM->Size(0)) {
    case 8:
      break;
    case 4:
      break;
    case 2:
      TestSVarExample1();
      break;
    case 1:
      TestSVarExample1();
      break;
    default:
      break;
    }
    break;
  case 8:
    switch (PPM->Size(0)) {
    case 8:
      break;
    case 4:
      TestSVarExample1();
      break;
    case 2:
      TestSVarExample1();
      break;
    case 1:
    default:
      break;
    }
    break;
  case 16:
    switch (PPM->Size(0)) {
    case 8:
      TestSVarExample1();
      break;
    case 4:
      TestSVarExample1();
      break;
    case 2:
      break;
    case 1:
    default:
      break;
    }
    break;
  case 1000:
    switch (PPM->Size(0)) {
    case 8:
      break;
    case 4:
      break;
    case 2:
      break;
    case 1:
      break;
    default:
      break;
    }
    break;
  }
}

class TestSLEstmNormalRVEvaluation : public TestSLEstimator {
protected:
  double TOL_MEAN = 4.0;

  double TOL_SVAR = 8.0;

  double TOL_SKEW = 16.0;

  double TOL_KURT = 64.0;
public:
  TestSLEstmNormalRVEvaluation() :
      TestSLEstimator(SLEstimatorConfig()
                          .WithDelayTotalUpdate(true)
                          .WithDelayParallelUpdate(true)
                          .WithInitLevel(2)
                          .WithMultiSampleFEMConfig(
                              MultiSampleFEMConfig("NormalLagrangeProtocol", "TODO")
                                  .WithPDESolverConfig(PDESolverConfig().WithDegree(1)))) {}
};

TEST_P(TestSLEstmNormalRVEvaluation, TestUpdate) {
  auto &aggregate = estimator.GetAggregate();
  for (auto &samples : {128, 256, 512, 1024, 2048, 4096, 8192}) {

    pout << estimator.GetAggregate() << endl;
    MemoryLogger::LogMemory();
    estimator.UpdateSampleAmount(samples);
    estimator.Method();
    estimator.GetMSFEM().GetMeshes().PrintInfo();
    pout << endl << estimator.GetAggregate() << endl;

    EXPECT_NEAR(0.0, aggregate.Q.GetComm().mean, TOL_MEAN * sqrt(1.0 / aggregate.ctr.Mcomm));
    EXPECT_NEAR(1.0, aggregate.Q.GetComm().sVar, TOL_SVAR * sqrt(1.0 / aggregate.ctr.Mcomm));
    EXPECT_NEAR(0.0, aggregate.Q.GetComm().skew, TOL_SKEW * sqrt(1.0 / aggregate.ctr.Mcomm));
    EXPECT_NEAR(3.0, aggregate.Q.GetComm().kurt, TOL_KURT * sqrt(1.0 / aggregate.ctr.Mcomm));
    EXPECT_DOUBLE_EQ(aggregate.Q.GetComm().mean, aggregate.Y.GetComm().mean);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetComm().sVar, aggregate.Y.GetComm().sVar);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetComm().skew, aggregate.Y.GetComm().skew);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetComm().kurt, aggregate.Y.GetComm().kurt);
#ifdef AGGREGATE_FOR_SOLUTION
    pout << "Mean Data after UpdateOnComm(): "
         << vec2str(aggregate.U.GetComm().mean.GetData().Data()) << endl
         << endl;
    for (row r = aggregate.U.GetComm().mean.rows(); r != aggregate.U.GetComm().mean.rows_end();
         r++) {
      EXPECT_NEAR(-r()[1], aggregate.U.GetComm().mean(r, 0),
                  TOL_MEAN * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(1.0, aggregate.U.GetComm().sVar(r, 0),
                  TOL_SVAR * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(0.0, aggregate.U.GetComm().skew(r, 0),
                  TOL_SKEW * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(3.0, aggregate.U.GetComm().kurt(r, 0),
                  TOL_KURT * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(-r()[1], aggregate.V.GetComm().mean(r, 0),
                  TOL_MEAN * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(1.0, aggregate.V.GetComm().sVar(r, 0),
                  TOL_SVAR * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(0.0, aggregate.V.GetComm().skew(r, 0),
                  TOL_SKEW * sqrt(1.0 / aggregate.ctr.Mcomm));
      EXPECT_NEAR(3.0, aggregate.V.GetComm().kurt(r, 0),
                  TOL_KURT * sqrt(1.0 / aggregate.ctr.Mcomm));
    }
#endif

    mout << "\t"
         << "calling UpdateParallel()" << endl
         << endl;
    estimator.UpdateParallel();
    pout << estimator.GetAggregate() << endl;

    EXPECT_NEAR(0.0, aggregate.Q.GetPara().mean, TOL_MEAN * sqrt(1.0 / aggregate.ctr.Mpara));
    EXPECT_NEAR(1.0, aggregate.Q.GetPara().sVar, TOL_SVAR * sqrt(1.0 / aggregate.ctr.Mpara));
    EXPECT_NEAR(0.0, aggregate.Q.GetPara().skew, TOL_SKEW * sqrt(1.0 / aggregate.ctr.Mpara));
    EXPECT_NEAR(3.0, aggregate.Q.GetPara().kurt, TOL_KURT * sqrt(1.0 / aggregate.ctr.Mpara));
    EXPECT_DOUBLE_EQ(aggregate.Q.GetPara().mean, aggregate.Y.GetPara().mean);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetPara().sVar, aggregate.Y.GetPara().sVar);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetPara().skew, aggregate.Y.GetPara().skew);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetPara().kurt, aggregate.Y.GetPara().kurt);
#ifdef AGGREGATE_FOR_SOLUTION
    pout << "Mean Data after UpdateParallel(): "
         << vec2str(aggregate.U.GetPara().mean.GetData().Data()) << endl
         << endl;
    for (row r = aggregate.U.GetPara().mean.rows(); r != aggregate.U.GetPara().mean.rows_end();
         r++) {
      EXPECT_NEAR(-r()[1], aggregate.U.GetPara().mean(r, 0),
                  TOL_MEAN * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(1.0, aggregate.U.GetPara().sVar(r, 0),
                  TOL_SVAR * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(0.0, aggregate.U.GetPara().skew(r, 0),
                  TOL_SKEW * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(3.0, aggregate.U.GetPara().kurt(r, 0),
                  TOL_KURT * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(-r()[1], aggregate.V.GetPara().mean(r, 0),
                  TOL_MEAN * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(1.0, aggregate.V.GetPara().sVar(r, 0),
                  TOL_SVAR * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(0.0, aggregate.V.GetPara().skew(r, 0),
                  TOL_SKEW * sqrt(1.0 / aggregate.ctr.Mpara));
      EXPECT_NEAR(3.0, aggregate.V.GetPara().kurt(r, 0),
                  TOL_KURT * sqrt(1.0 / aggregate.ctr.Mpara));
    }
#endif

    mout << "\t"
         << "calling UpdateTotal()" << endl
         << endl;
    estimator.UpdateTotal();
    pout << estimator.GetAggregate() << endl;

    EXPECT_NEAR(0.0, aggregate.Q.GetTotal().mean, TOL_MEAN * sqrt(1.0 / aggregate.ctr.M));
    EXPECT_NEAR(1.0, aggregate.Q.GetTotal().sVar, TOL_SVAR * sqrt(1.0 / aggregate.ctr.M));
    EXPECT_NEAR(0.0, aggregate.Q.GetTotal().skew, TOL_SKEW * sqrt(1.0 / aggregate.ctr.M));
    EXPECT_NEAR(3.0, aggregate.Q.GetTotal().kurt, TOL_KURT * sqrt(1.0 / aggregate.ctr.M));
    EXPECT_DOUBLE_EQ(aggregate.Q.GetTotal().mean, aggregate.Y.GetTotal().mean);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetTotal().sVar, aggregate.Y.GetTotal().sVar);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetTotal().skew, aggregate.Y.GetTotal().skew);
    EXPECT_DOUBLE_EQ(aggregate.Q.GetTotal().kurt, aggregate.Y.GetTotal().kurt);
#ifdef AGGREGATE_FOR_SOLUTION
    pout << "Mean Data after UpdateTotal(): "
         << vec2str(aggregate.U.GetTotal().mean.GetData().Data()) << endl
         << endl;
    for (row r = aggregate.U.GetTotal().mean.rows(); r != aggregate.U.GetTotal().mean.rows_end();
         r++) {
      EXPECT_NEAR(-r()[1], aggregate.U.GetTotal().mean(r, 0),
                  TOL_MEAN * sqrt(1.0 / aggregate.ctr.M));
      EXPECT_NEAR(1.0, aggregate.U.GetTotal().sVar(r, 0), TOL_SVAR * sqrt(1.0 / aggregate.ctr.M));
      EXPECT_NEAR(0.0, aggregate.U.GetTotal().skew(r, 0), TOL_SKEW * sqrt(1.0 / aggregate.ctr.M));
      //      EXPECT_NEAR(3.0, aggregate.U.GetTotal().kurt(r, 0), TOL_KURT * sqrt(1.0 /
      //      aggregate.ctr.M));
      EXPECT_NEAR(-r()[1], aggregate.V.GetTotal().mean(r, 0),
                  TOL_MEAN * sqrt(1.0 / aggregate.ctr.M));
      EXPECT_NEAR(1.0, aggregate.V.GetTotal().sVar(r, 0), TOL_SVAR * sqrt(1.0 / aggregate.ctr.M));
      EXPECT_NEAR(0.0, aggregate.V.GetTotal().skew(r, 0), TOL_SKEW * sqrt(1.0 / aggregate.ctr.M));
      //      EXPECT_NEAR(3.0, aggregate.V.GetTotal().kurt(r, 0), TOL_KURT * sqrt(1.0 /
      //      aggregate.ctr.M));
    }
#endif
  }
}

INSTANTIATE_TEST_SUITE_P(TestSLEstimator, TestSLEstmProcEvaluation, Values(4, 8, 16, 32));

INSTANTIATE_TEST_SUITE_P(TestSLEstimator, TestSLEstmModuloIndexEvaluation, Values(4, 8, 16, 32));

INSTANTIATE_TEST_SUITE_P(TestSLEstimator, TestSLEstmZeroEvaluation, Values(4, 8, 16, 32));

INSTANTIATE_TEST_SUITE_P(TestSLEstimator, TestSLEstmNormalRVEvaluation, Values(0));

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithParallelListeners()
                     .WithConfigEntry("Distribution", "Stripes_y")
                     .WithConfigEntry("SampleCounterVerbose", 1)
                     .WithConfigEntry("SLEstimatorVerbose", 1)
                     .WithScreenLogging()
                     .WithRandomInitialized()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}