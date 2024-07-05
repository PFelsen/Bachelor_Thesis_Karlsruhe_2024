#include "TestDistribution.hpp"

TEST_P(TestDistribution, TestCellCount) {
  for (int commSplit = 0; commSplit <= maxCommSplit; commSplit++) {
    int cells1 = (*meshes)[{pLevel, -1, 0, commSplit}].CellCount();
    int cells2 = ceil((double)(*meshes)[{pLevel, -1, 0, commSplit}].CellCountGeometry()
                      / PPM->Size(commSplit));
    EXPECT_EQ(cells1, cells2);
    cells1 = (*meshes)[{level, -1, 0, commSplit}].CellCount();
    cells2 = ceil((double)(*meshes)[{level, -1, 0, commSplit}].CellCountGeometry()
                  / PPM->Size(commSplit));
    EXPECT_EQ(cells1, cells2);
    vout(1) << "level=" << level << " commSplit=" << commSplit << endl;
    vout(1) << cells1 << " " << cells2 << endl;
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("TestDistributionVerbose", 1)
                     .WithConfigEntry("DistributionVerbose", 0)
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}