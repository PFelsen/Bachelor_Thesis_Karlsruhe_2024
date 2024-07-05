#include "TestTimeIntegrator.hpp"
#include "GMRES.hpp"
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

  std::vector<TestAssemble *> assembles{
      new TestPAssemble<0>(), new TestPAssemble<1>(), new TestPAssemble<2>(),
      new TestPAssemble<3>(), new TestPAssemble<4>(), new TestPAssemble<5>(),
  };

  std::shared_ptr<const LagrangeDiscretization> disc;

  Vector uInt;

  Vector uEx;

  TestTimeIntegrator() :
      meshes(MeshesCreator("Interval").WithPLevel(3).WithLevel(3).CreateUnique()),
      timeInt(TimeIntegratorCreator(GetParam().solverID)
                  .WithLinearSolver(new GMRES(GetPC("SuperLU")))
                  .WithLinearSolver(new GMRES(GetPC("SuperLU")))
                  .WithLinearSolver(new GMRES(GetPC("SuperLU")))
                  .WithLinearSolver(new GMRES(GetPC("SuperLU")))
                  .CreateUnique()),
      disc(std::make_shared<const LagrangeDiscretization>(*meshes, 0)), uInt(disc), uEx(disc) {
    uInt = 1.0;
    uInt.MakeAdditive();
    uEx = 1.0;
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
                                TestParams{DI_RK_ORDER3, 3} //    TestParams{ARNOLDI2, 2},
                                                            //    TestParams{EXPONENTIAL, 0},
                                ));

TEST_P(TestTimeIntegrator, TestGetOrder) { EXPECT_EQ(timeInt->GetOrder(), GetParam().solverOder); }

TEST_P(TestTimeIntegrator, TestIntegratePolynomial) {
  for (int degree = 1; degree <= timeInt->GetOrder(); degree++) {
    timeInt->PrintInfo();
    timeInt->Method(assembles[degree], uInt);
    assembles[degree]->Solution(assembles[degree]->LastTStep(), uEx);
    EXPECT_DOUBLE_EQ(uInt.norm(), uEx.norm());
  }

  for (int degree = timeInt->GetOrder() + 1; degree <= 5; degree++) {
    timeInt->Method(assembles[degree], uInt);
    assembles[degree]->Solution(assembles[degree]->LastTStep(), uEx);
    if (degree != 4)
      EXPECT_NE(uInt.norm(), uEx.norm()); // have to exclude tihs specific order 3 bc.
    // the nodes/weights are gaußian and therefore for a problem that only depends on time
    // integrates exactly
    //    else Warning("degree " + to_string(degree) + " excluded from test")
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("TimeIntegratorVerbose", 3)
                     .WithConfigEntry("TimeSeriesVerbose", 1)
                     .WithConfigEntry("LinearVerbose", 0)
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
