#include "TestCell.hpp"
#include "TestEnvironment.hpp"

class CellInterfaceTest : public CellTest {
public:
  CellInterfaceTest() : CellTest(getReferenceParameters(GetParam().celltype)) {}
};

// TODO: 07.09.21: Methods now access vector via [] and won't throw. How to test invalid values?

/*
TEST_P(CellInterfaceTest, TestNonexistingCorners) {
    EXPECT_ANY_THROW(c->Corner(GetParam().corners.size()));
}


TEST_P(CellInterfaceTest, TestNonexistingEdges) {
    EXPECT_ANY_THROW(c->Edge(GetParam().edges.size()));
}

TEST_P(CellInterfaceTest, TestNonexistinFaces) {
    EXPECT_ANY_THROW(c->z(Face(GetParam().faces.size())));
}*/

// INSTANTIATE_TEST_SUITE_P(InterfaceTests, CellInterfaceTest, ValuesIn(referenceCells));

int main(int argc, char **argv) {
  constructTestData(cellTestData);
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}