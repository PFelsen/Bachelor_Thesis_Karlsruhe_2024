#include "gtest/gtest.h"

#include "TestBasicAlgebra.hpp"
#include "basicalgebra/CTensor.hpp"

class CTensorTest : public BasicAlgebraTensorTest {
protected:
  CTensor A, B;
  RTensor A_d;
  std::complex<double> value;
  double value_d;
  int value_i;
public:
  CTensorTest() : BasicAlgebraTensorTest() {
    mpp_ba::SetTolerance(5e-14);
    A.resize(FirstComponents, SecondComponents, ThirdComponents);
    B.resize(FirstComponents, SecondComponents, ThirdComponents);
    A_d.resize(FirstComponents, SecondComponents, ThirdComponents);

    // A(RandomDouble(), FirstComponents, SecondComponents, ThirdComponents);

    for (int i = 0; i < FirstComponents; ++i) {
      for (int j = 0; j < SecondComponents; ++j) {
        for (int k = 0; k < ThirdComponents; ++k) {
          B(i, j, k) = RandomComplex();
          A_d(i, j, k) = RandomDouble();
        }
      }
    }

    value = RandomNonZeroComplex();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();
  }
};

TEST_P(CTensorTest, ConstructorTest) {
  CTensor W0(FirstComponents, SecondComponents, ThirdComponents);
  CTensor W1(value, FirstComponents, SecondComponents, ThirdComponents);
  CTensor W2(value_d, FirstComponents, SecondComponents, ThirdComponents);
  CTensor W3(value_i, FirstComponents, SecondComponents, ThirdComponents);
  CTensor W4(B);
  CTensor W5(A_d);

  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(W0(i, j, k), std::complex<double>(0.0));
        EXPECT_COMPLEX_EQ(W1(i, j, k), value);
        EXPECT_COMPLEX_EQ(W2(i, j, k), std::complex<double>(value_d));
        EXPECT_COMPLEX_EQ(W3(i, j, k), std::complex<double>(value_i));
        EXPECT_COMPLEX_EQ(W4(i, j, k), B(i, j, k));
        EXPECT_COMPLEX_EQ(W5(i, j, k), std::complex<double>(A_d(i, j, k)));
      }
    }
  }
}

TEST_P(CTensorTest, SummationTest) {
  CTensor AB = A + B;
  CTensor W1 = A + A_d;
  CTensor W2 = A_d + B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(AB(i, j, k), A(i, j, k) + B(i, j, k));
        EXPECT_COMPLEX_EQ(W1(i, j, k), A(i, j, k) + std::complex<double>(A_d(i, j, k)));
        EXPECT_COMPLEX_EQ(W2(i, j, k), std::complex<double>(A_d(i, j, k)) + B(i, j, k));
      }
    }
  }
}

TEST_P(CTensorTest, DifferenceTest) {
  CTensor AB = A - B;
  CTensor W1 = A - A_d;
  CTensor W2 = A_d - B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(AB(i, j, k), A(i, j, k) - B(i, j, k));
        EXPECT_COMPLEX_EQ(W1(i, j, k), A(i, j, k) - std::complex<double>(A_d(i, j, k)));
        EXPECT_COMPLEX_EQ(W2(i, j, k), std::complex<double>(A_d(i, j, k)) - B(i, j, k));
      }
    }
  }
}

TEST_P(CTensorTest, AdditiveInversTest) {
  CTensor W = -B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(W(i, j, k), -B(i, j, k));
      }
    }
  }
}

TEST_P(CTensorTest, ScalarMultiplicationTest) {
  CTensor W1 = value * B;
  CTensor W2 = B * value;
  CTensor W3 = value_d * B;
  CTensor W4 = B * value_d;
  CTensor W5 = value_i * B;
  CTensor W6 = B * value_i;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(W1(i, j, k), value * B(i, j, k));
        EXPECT_COMPLEX_EQ(W2(i, j, k), value * B(i, j, k));
        EXPECT_COMPLEX_EQ(W3(i, j, k), value_d * B(i, j, k));
        EXPECT_COMPLEX_EQ(W4(i, j, k), value_d * B(i, j, k));
        EXPECT_COMPLEX_EQ(W5(i, j, k), value_i * B(i, j, k));
        EXPECT_COMPLEX_EQ(W6(i, j, k), value_i * B(i, j, k));
      }
    }
  }
}

TEST_P(CTensorTest, ScalarDivisionTest) {
  CTensor W1 = B / value;
  CTensor W2 = B / value_d;
  CTensor W3 = B / value_i;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(W1(i, j, k), B(i, j, k) / value);
        EXPECT_COMPLEX_EQ(W2(i, j, k), B(i, j, k) / value_d);
        EXPECT_COMPLEX_EQ(W3(i, j, k), B(i, j, k) / value_i);
      }
    }
  }
}

TEST_P(CTensorTest, TensorMultiplicationTest) {
  CTensor W0(A, B);
  CTensor W1(B, A);
  CTensor W2(B, B);
  CTensor W3, W4, W5, W6, W7;
  W3 = A_d * B;
  W4 = B * A_d;
  W5 = A * B;
  W6 = B * A;
  W7 = B * B;
  for (int i = 0; i < FirstComponents; ++i) {
    for (int j = 0; j < SecondComponents; ++j) {
      for (int k = 0; k < ThirdComponents; ++k) {
        EXPECT_COMPLEX_EQ(W0(i, j, k), A(i, j, k) * B(i, j, k));
        EXPECT_COMPLEX_EQ(W1(i, j, k), A(i, j, k) * B(i, j, k));
        EXPECT_COMPLEX_EQ(W2(i, j, k), B(i, j, k) * B(i, j, k));

        EXPECT_COMPLEX_EQ(W3(i, j, k), std::complex<double>(A_d(i, j, k)) * B(i, j, k));
        EXPECT_COMPLEX_EQ(W4(i, j, k), std::complex<double>(A_d(i, j, k)) * B(i, j, k));
        EXPECT_COMPLEX_EQ(W5(i, j, k), A(i, j, k) * B(i, j, k));
        EXPECT_COMPLEX_EQ(W6(i, j, k), A(i, j, k) * B(i, j, k));
        EXPECT_COMPLEX_EQ(W7(i, j, k), B(i, j, k) * B(i, j, k));
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    BasicAlgebraTensorTest, CTensorTest,
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