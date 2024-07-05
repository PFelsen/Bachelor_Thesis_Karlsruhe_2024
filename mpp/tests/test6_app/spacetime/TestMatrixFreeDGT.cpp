#include <map>
#include <vector>
#include "DGViscoAcousticDGDG_GLGL.hpp"
#include "Results.hpp"
#include "SpaceTime.hpp"
#include "SpaceTimeDiscretizations.hpp"
#include "TestEnvironment.hpp"

const double test_tol = 1e-10;

class MatrixFreeDGT : public Test {
protected:
  std::unique_ptr<Meshes> STM;
  ProblemBase *problem;

  void SetUp() override {
    STM = MeshesCreator().WithMeshName("ST_squ_simple").WithPLevel(1).WithLevel(1).CreateUnique();
  }

  void TearDown() override {}

  /*
    template<typename SystemOperator>
    void testMatrixFree(std::shared_ptr<STDiscretization> &disc){
      SpaceTimeMatrixGraph cellmg(*mesh, *disc, false);
      auto mg = VectorMatrixBase(cellmg);
      auto assemble = TDGViscoAcousticAssemble_DGT{*disc, problem};
      Preconditioner *PC = GetPC("PointBlockJacobi");
      Vector RHS(0.0, mg);
      Vector v(0.0, mg);
      fillSpaceTimeVector(v, problem.problemDim(),
                          *disc,
                          [](const Point& p, const tcell& c, int i){
        return p[i];
      });
      Vector B1(0.0,mg);
      Vector B2(0.0,mg);
      Vector B2_s(0.0,mg);
      Matrix B = Matrix(mg);
      B = 0;
      {

        mout.StartBlock("FullAssemble_OFFLINE");
        assemble.System(B, RHS);
        PC->Construct(B);
        mout << "Memory: " << B.Size()* 8.0 / 1024.0 / 1024.0 << " MB" << endl;
        mout.EndBlock();
        mout.StartBlock("FullAssemble__ONLINE");
        B1 = B*v;
        mout << "|Ax| = " << B1.norm() << endl;
        B2 = (*PC) * B1;
        B2_s = B2;
        mout << "|BA1| =  " << B2.norm() << endl;
        B1 = B * B2;
        mout << "|ABA1| =  " << B1.norm() << endl;
        mout.EndBlock();
        mout << endl;
      }

      mout.StartBlock("PartAssemble_OFFLINE");
      auto *B_mf = new SystemOperator(assemble, mg, true);
      mout << B_mf->Info() << endl;
      mout.EndBlock();
      mout.StartBlock("PartAssemble__ONLINE");
      Vector B1_mfree = (*B_mf) * v;
      mout << "|Ax| = " << B1_mfree.norm() << endl;
      B2 = (*PC) * B1_mfree;
      mout << "|BA1| =  " << B2.norm() << endl;
      B1_mfree = (*B_mf) * B2;
      double dd = (B2-B2_s).norm();
      mout << "|BA1 - BA1|.norm(): " << dd << endl;
      mout << "|ABA1| =  " << B1_mfree.norm() << endl;
      double diff = (B1-B1_mfree).norm();
      mout << "|ABA1 - ABA1|.norm(): " << diff << endl;
      mout.EndBlock();
      mout << endl;

      mout.StartBlock("NoAssemble_OFFLINE");
      B_mf = new SystemOperator(assemble, mg, false);
      mout << B_mf->Info() << endl;
      mout.EndBlock();
      mout.StartBlock("NoAssemble__ONLINE");
      B1_mfree = (*B_mf) * v;
      mout << "|Ax| = " << B1_mfree.norm() << endl;
      B2 = (*PC) * B1_mfree;
      mout << "|BA1| =  " << B2.norm() << endl;
      dd = (B2-B2_s).norm();
      mout << "|BA1 - BA1|.norm(): " << dd << endl;
      B1_mfree = (*B_mf) * B2;
      mout << "|ABA1| =  " << B1_mfree.norm() << endl;
      mout << B1-B1_mfree << endl;
      diff = (B2-B2_s).norm();
      mout << "|ABA1 - ABA1|.norm(): " << diff << endl;
      mout.EndBlock();

    }*/
};

/*
TEST_F(MatrixFreeDGT, MatrixFree_EQ) {

  auto disc = createDiscretization(1,1,problem.problemDim(),
                                   *mesh,true, "EQ");

  testMatrixFree<TDGViscoAcoustic_DGT_SystemOperator>(disc);
}*/

TEST_F(MatrixFreeDGT, MatrixFree_GL) {

  auto disc = CreateSTDiscretizationShared(*STM, 1, 1, problem->Dim(), true, "GLGL");

  // testMatrixFree<TDGViscoAcoustic_DGT_SystemOperator_GLGL>(disc);
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithScreenLogging().WithPPM().WithoutDefaultConfig();
  std::map<string, string> map{{"Distribution", "deformed_optimized"},
                               {"Verbose", "5"},
                               {"ConfigVerbose", "0"}};

  Config::Initialize(map);

  return mppTest.RUN_ALL_MPP_TESTS();
}
