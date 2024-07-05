#include <memory>

#include "TestEnvironment.hpp"
#include "spectrum/EigenvalueExample.hpp"

const double ev_tol = 1e-3;

class EigenvalueExampleTest : public ::testing::Test {
protected:
  int N_ev = 6;
  double shift = 2.0;
  double rho = 127.0;

  std::unique_ptr<Meshes> meshes = nullptr;
  std::unique_ptr<EigenvalueExampleAssemble<true>> assemble = nullptr;
  std::shared_ptr<IDiscretization> disc_ev = nullptr;
  std::shared_ptr<IAIDiscretization> ia_disc_ev = nullptr;

  std::vector<IAInterval> eigenvalues{};

  EigenvalueExampleTest() {
    meshes = MeshesCreator("Square2Triangles").WithPLevel(1).WithLevel(5).CreateUnique();
    meshes->PrintInfo();

    assemble = std::make_unique<EigenvalueExampleAssemble<true>>(shift);
    disc_ev = assemble->CreateEVDisc(*meshes);
    ia_disc_ev = assemble->CreateIAEVDisc(*meshes);

    eigenvalues = ComputeEigenvaluesBaseProblem(N_ev);
  }
};

TEST_F(EigenvalueExampleTest, RayleighRitzTest) {
  IAUpperEigenvalueBounds Lambda =
      ApplyRayleighRitz<true>(*meshes, *assemble, disc_ev, ia_disc_ev, N_ev, shift);
  for (int i = 0; i < Lambda.size(); ++i)
    EXPECT_LE(sup(eigenvalues[i]), Lambda[i]);
}

TEST_F(EigenvalueExampleTest, LehmannGoerischTest) {
  IALowerEigenvalueBounds lambda =
      ApplyLehmannGoerisch<true>(*meshes, *assemble, disc_ev, ia_disc_ev, N_ev, shift, rho);
  for (int i = 0; i < lambda.size(); ++i)
    EXPECT_LE(lambda[i], inf(eigenvalues[i]));
}

TEST_F(EigenvalueExampleTest, HomotopyTest) {
  double final_step =
      ApplyHomotopyMethod<true>(*meshes, *assemble, disc_ev, ia_disc_ev, N_ev, shift, rho);
  IAEigenvalueEnclosures IALambda = ComputeEnclosures<true>(*meshes, *assemble, disc_ev, ia_disc_ev,
                                                            N_ev, shift, rho, final_step);
  for (int i = 0; i < IALambda.size(); ++i) {
    EXPECT_LE(IALambda.getLowerBound(i), IALambda.getUpperBound(i));
    EXPECT_NEAR(IALambda.getLowerBound(i), IALambda.getUpperBound(i), ev_tol);
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfPath(std::string(ProjectMppDir) + "/conf/")
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}
