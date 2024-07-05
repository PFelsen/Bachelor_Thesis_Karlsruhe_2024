#include <LagrangeDiscretization.hpp>
#include <MeshesCreator.hpp>

#include "LagrangeTransfer.hpp"
#include "TestTransfers.hpp"

class TestLagrangeTransfer : public TestWithParam<TransferParam> {
protected:
  std::unique_ptr<Meshes> meshes;
  std::shared_ptr<const LagrangeDiscretization> coarseDisc;
  std::shared_ptr<const LagrangeDiscretization> fineDisc;
  std::unique_ptr<Vector> coarse;
  std::unique_ptr<Vector> fine;
  std::unique_ptr<ITransfer> T;

  TestLagrangeTransfer() {
    meshes = std::unique_ptr<Meshes>(
        MeshesCreator(GetParam().meshName).WithLevel(GetParam().fineLevel).Create());

    coarseDisc = std::make_shared<LagrangeDiscretization>(*meshes, GetParam().coarseDegree, 1);
    coarse = std::make_unique<Vector>(coarseDisc, GetParam().coarseLevel);
    fineDisc = std::make_shared<LagrangeDiscretization>(*meshes, GetParam().fineDegree, 1);
    fine = std::make_unique<Vector>(fineDisc, GetParam().fineLevel);

    T = std::make_unique<ITransfer>(new LagrangeTransfer(*coarse, *fine));
  }

  void SetUp() override { setSeed((unsigned)std::time(nullptr)); }
};

TEST_P(TestLagrangeTransfer, ConstantProject) {
  double d = RandomDouble();
  FillVectorConstant(*coarse, 0.0);
  FillVectorConstant(*fine, d);

  T->Project(*coarse, *fine);

  ASSERT_CONSTANTVECTOR(*coarse, d);
}

TEST_P(TestLagrangeTransfer, LinearProject) {
  double dmin = RandomDouble();
  double dmax = RandomDouble();
  FillVectorConstant(*coarse, 0.0);
  FillVectorLinear(*fine, dmin, dmax);

  T->Project(*coarse, *fine);

  ASSERT_LINEARVECTOR(*coarse, dmin, dmax);
}

TEST_P(TestLagrangeTransfer, ConstantRestrict) {
  double d = RandomDouble();
  FillVectorConstant(*coarse, 0.0);
  FillVectorConstant(*fine, d);

  T->Restrict(*coarse, *fine);

  ASSERT_RESTRICTION(*coarse, *fine, T->AsMatrix());
}

TEST_P(TestLagrangeTransfer, LinearRestrict) {
  if (PPM->size() > 0) return;
  double dmin = RandomDouble();
  double dmax = RandomDouble();
  FillVectorConstant(*coarse, 0.0);
  FillVectorLinear(*fine, dmin, dmax);

  T->Restrict(*coarse, *fine);

  ASSERT_RESTRICTION(*coarse, *fine, T->AsMatrix());
}

TEST_P(TestLagrangeTransfer, ConstantElongate) {
  double d = RandomDouble();
  FillVectorConstant(*coarse, d);
  FillVectorConstant(*fine, 0.0);

  T->Prolongate(*coarse, *fine);

  ASSERT_CONSTANTVECTOR(*fine, d);
}

TEST_P(TestLagrangeTransfer, LinearElongate) {
  double dmin = RandomDouble();
  double dmax = RandomDouble();
  FillVectorLinear(*coarse, dmin, dmax);
  FillVectorConstant(*fine, 0.0);

  T->Prolongate(*coarse, *fine);

  ASSERT_LINEARVECTOR(*fine, dmin, dmax);
}

INSTANTIATE_TEST_SUITE_P(Interval, TestLagrangeTransfer, ValuesIn(TransferParameters("Interval")));
#if SpaceDimension >= 2
INSTANTIATE_TEST_SUITE_P(Triangle, TestLagrangeTransfer, ValuesIn(TransferParameters("Triangle")));
INSTANTIATE_TEST_SUITE_P(Square, TestLagrangeTransfer, ValuesIn(TransferParameters("Square")));
INSTANTIATE_TEST_SUITE_P(FourTriangles, TestLagrangeTransfer,
                         ValuesIn(TransferParameters("Square4Triangles")));
#endif
#if SpaceDimension >= 3
INSTANTIATE_TEST_SUITE_P(Tetrahedron, TestLagrangeTransfer,
                         ValuesIn(TransferParameters("Tetrahedron")));
INSTANTIATE_TEST_SUITE_P(Hexahedron, TestLagrangeTransfer,
                         ValuesIn(TransferParameters("Hexahedron")));
#endif


int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}