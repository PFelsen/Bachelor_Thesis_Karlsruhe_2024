//
// Created by cr on 14.11.21.
//

#include "TestTimeIntegratorEOC.h"
#include <DGDiscretization.hpp>
#include <GMRES.hpp>
#include <MeshesCreator.hpp>
#include <TimeIntegrator.hpp>
#include "TestEnvironment.hpp"

struct TestParams {
  TIMEINTEGRATOR solverID;
  int solverOder;
};

std::ostream &operator<<(std::ostream &s, const TestParams &testParams) {
  return s << "{Solver ID: " << testParams.solverID << ", Solver Order: " << testParams.solverOder
           << "}" << endl;
}

class TestTimeIntegrator : public TestWithParam<TestParams> {
protected:
  std::unique_ptr<Meshes> meshes;

  std::unique_ptr<TimeIntegrator> timeInt;

  std::shared_ptr<const DGDiscretization> disc;

  std::vector<TestAssembleEOC *> assembles{
      new TestAssembleEOC(1), new TestAssembleEOC(2), new TestAssembleEOC(3),
      new TestAssembleEOC(4), new TestAssembleEOC(5), new TestAssembleEOC(6),
  };
  Vector uInt;

  Vector uEx;

  TestTimeIntegrator() :
      meshes(MeshesCreator("Interval").WithPLevel(4).WithLevel(4).CreateUnique()),
      timeInt(TimeIntegratorCreator(GetParam().solverID)
                  .WithLinearSolver(new GMRES(GetPC("PointBlockJacobi")))
                  .WithLinearSolver(new GMRES(GetPC("PointBlockJacobi")))
                  .WithLinearSolver(new GMRES(GetPC("PointBlockJacobi")))
                  .WithLinearSolver(new GMRES(GetPC("PointBlockJacobi")))
                  .CreateUnique()),
      disc(std::make_shared<const DGDiscretization>(*meshes, 0, 2)), uInt(disc), uEx(disc) {
    uInt = 1.0;
    uInt.MakeAdditive();
    uEx = 1.0;
    uEx.MakeAdditive();
  }

  void TearDown() override {
    for (auto &assemble : assembles)
      delete assemble;
    assembles.clear();
  }
};

INSTANTIATE_TEST_SUITE_P(TestTimeIntegrator, TestTimeIntegrator,
                         Values(TestParams{EXPLICIT_EULER, 1}, TestParams{HEUN, 2},
                                TestParams{RUNGE, 3}, TestParams{RUNGE_KUTTA, 4},
                                TestParams{IMPLICIT_EULER, 1}, TestParams{IMPLICIT_MIDPOINT, 2},
                                TestParams{CRANK_NICOLSON, 2}, TestParams{DI_RK_ORDER2, 2},
                                TestParams{DI_RK_ORDER3, 3}, TestParams{EXPONENTIAL_EULER, 1},
                                TestParams{EXPONENTIAL_MIDPOINT, 2},
                                TestParams{EXPONENTIAL_SHIFT_EULER, 1},
                                TestParams{EXPONENTIAL_SHIFT_MIDPOINT, 2} //
                                //        TestParams{EXPONENTIAL_EXTENDED_EULER, 1}//

                                //    TestParams{ARNOLDI2, 2},
                                //    TestParams{EXPONENTIAL, 0},
                                ));

TEST_P(TestTimeIntegrator, TestGetOrder) { EXPECT_EQ(timeInt->GetOrder(), GetParam().solverOder); }

TEST_P(TestTimeIntegrator, TestIntegratePolynomial) {
  timeInt->PrintInfo();
  std::vector<Vector> uInts;
  for (int i = 0; i < assembles.size(); i++) {
    timeInt->Method(assembles[i], uInt);
    uInts.push_back(uInt);
  }
  std::string s;
  for (int i = 0; i < assembles.size(); i++) {
    assembles[i]->Solution(assembles[i]->LastTStep(), uEx);
    //        uInts[i].MakeAdditive();
    s += std::to_string((uEx - uInts[i]).norm() / uEx.norm()) + " ";
  }
  mout << s << endl;
  std::vector<double> diffNorms(assembles.size() - 1, 0);
  for (int i = 1; i < assembles.size(); i++) {
    diffNorms[i - 1] = (uInts[i] - uInts[i - 1]).norm();
  }
  double sum = 0;
  for (int i = 0; i < assembles.size() - 2; ++i) {
    //    MOUT(log(diffNorms[i] / diffNorms[i + 1]) / log(2));
    sum += log(diffNorms[i] / diffNorms[i + 1]) / log(2);
  }
  MOUT(sum / (assembles.size() - 2))
  EXPECT_NEAR(timeInt->GetOrder(), sum / (assembles.size() - 2), timeInt->GetOrder() * 0.1);
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("TimeIntegratorVerbose", -3)
                     .WithConfigEntry("TimeSeriesVerbose", 1)
                     .WithConfigEntry("LinearVerbose", -1)
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}