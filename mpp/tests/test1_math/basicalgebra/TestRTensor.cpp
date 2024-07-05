#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/RTensor.hpp"

class RTensorTest : public BasicAlgebraTensorTest {
protected:
  RTensor A, B;
  double value;
  int value_i;
public:
  RTensorTest() : BasicAlgebraTensorTest() {
    A.resize(FirstComponents, SecondComponents, ThirdComponents);
    B.resize(FirstComponents, SecondComponents, ThirdComponents);

    // A(RandomDouble(), FirstComponents, SecondComponents, ThirdComponents);

    for (int i = 0; i < FirstComponents; ++i) {
      for (int j = 0; j < SecondComponents; ++j) {
        for (int k = 0; k < ThirdComponents; ++k) {
          B(i, j, k) = RandomDouble();
        }
      }
    }

    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(RTensorTest, ConstructorTest) {
  RTensor W0(FirstComponents, SecondComponents, ThirdComponents);
  RTensor W1(value, FirstComponents, SecondComponents, ThirdComponents);
  RTensor W2(value_i, FirstComponents, SecondComponents, ThirdComponents);
  ;
  RTensor W3(B);

  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(W0(i, j, k), 0.0);
        EXPECT_DOUBLE_EQ(W1(i, j, k), value);
        EXPECT_DOUBLE_EQ(W2(i, j, k), value_i);
        EXPECT_DOUBLE_EQ(W3(i, j, k), B(i, j, k));
      }
    }
  }
}

TEST_P(RTensorTest, SummationTest) {
  RTensor AB = A + B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(AB(i, j, k), A(i, j, k) + B(i, j, k));
      }
    }
  }
}

TEST_P(RTensorTest, DifferenceTest) {
  RTensor AB = A - B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(AB(i, j, k), A(i, j, k) - B(i, j, k));
      }
    }
  }
}

TEST_P(RTensorTest, AdditiveInversTest) {
  RTensor W = -B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(W(i, j, k), -B(i, j, k));
      }
    }
  }
}

TEST_P(RTensorTest, ScalarMultiplicationTest) {
  RTensor W1 = value * B;
  RTensor W2 = B * value;
  RTensor W3 = value_i * B;
  RTensor W4 = B * value_i;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(W1(i, j, k), value * B(i, j, k));
        EXPECT_DOUBLE_EQ(W2(i, j, k), value * B(i, j, k));
        EXPECT_DOUBLE_EQ(W3(i, j, k), value_i * B(i, j, k));
        EXPECT_DOUBLE_EQ(W4(i, j, k), value_i * B(i, j, k));
      }
    }
  }
}

TEST_P(RTensorTest, ScalarDivisionTest) {
  RTensor W1 = B / value;
  RTensor W2 = B / value_i;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(W1(i, j, k), B(i, j, k) / value);
        EXPECT_DOUBLE_EQ(W2(i, j, k), B(i, j, k) / value_i);
      }
    }
  }
}

TEST_P(RTensorTest, TensorMultiplicationTest) {
  RTensor W0(A, B);
  RTensor W1(B, A);
  RTensor W2(B, B);
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_DOUBLE_EQ(W0(i, j, k), A(i, j, k) * B(i, j, k));
        EXPECT_DOUBLE_EQ(W1(i, j, k), A(i, j, k) * B(i, j, k));
        EXPECT_DOUBLE_EQ(W2(i, j, k), B(i, j, k) * B(i, j, k));
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    BasicAlgebraTensorTest, RTensorTest,
    Values(BasicAlgebraTensorTestParameter{1, 1, 1}, BasicAlgebraTensorTestParameter{2, 1, 1},
           BasicAlgebraTensorTestParameter{3, 1, 1},

           BasicAlgebraTensorTestParameter{2, 1, 1}, BasicAlgebraTensorTestParameter{2, 2, 1},
           BasicAlgebraTensorTestParameter{2, 3, 1},

           BasicAlgebraTensorTestParameter{3, 3, 1}, BasicAlgebraTensorTestParameter{3, 3, 2},
           BasicAlgebraTensorTestParameter{3, 3, 3}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}