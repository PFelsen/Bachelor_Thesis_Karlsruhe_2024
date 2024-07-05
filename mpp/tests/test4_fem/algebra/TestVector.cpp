#include "TestVector.hpp"

TEST_F(TestVector, TestDoubleAccumulationOverlap) {
  vector = 1.0;
  vector.Accumulate();
  testVectorOverlap(vector, 1.0);
}

TEST_F(TestVector, TestDoubleAccumulationNonOverlap) {
  vector = 1.0;
  vector.Accumulate();
  testVectorNonOverlap(vector, 1.0);
}

TEST_F(TestVector, TestAddingWithoutOverlap) {
  vector = 1.0;
  vector.MakeAdditive();
  Vector vector2(vector);
  vector2 = 4.0;
  testAdditionWithoutOverlap(vector, vector2);
}

TEST_F(TestVector, TestNormWithoutOverlap) {
  vector = 1.0;
  testNormWithoutOverlap(vector);
}

TEST_F(TestVector, TestHadamard) {
  Vector product(0.0, vector);
  product = HadamardProductVector(vector, product);
  for (row r = product.rows(); r != product.rows_end(); ++r) {
    EXPECT_EQ(product(r, 0), 0.0);
  }
  product = HadamardProductVector(vector, vector);
  for (row r = product.rows(); r != product.rows_end(); ++r) {
    EXPECT_EQ(product(r, 0), 1.0);
  }
}

TEST_F(TestVector, TestSqrtComponents) {
  Vector root(4.0, vector);
  root = ComponentSqrt(root);
  for (row r = root.rows(); r != root.rows_end(); ++r) {
    EXPECT_EQ(root(r, 0), 2.0);
  }
  root = ComponentSqrt(vector);
  for (row r = root.rows(); r != root.rows_end(); ++r) {
    EXPECT_EQ(root(r, 0), 1.0);
  }
}

TEST_F(TestVector, TestComponentDivide) {
  Vector division(0.5, vector);
  division = ComponentDivide(vector, division);
  for (row r = division.rows(); r != division.rows_end(); ++r) {
    EXPECT_EQ(division(r, 0), 2.0);
  }
  division = ComponentDivide(vector, vector);
  for (row r = division.rows(); r != division.rows_end(); ++r) {
    EXPECT_EQ(division(r, 0), 1.0);
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("DistributionVerbose", 0)
                     .WithConfigEntry("MeshesVerbose", 1)
                     .WithConfigEntry("MeshVerbose", 2)
                     .WithParallelListeners()
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}