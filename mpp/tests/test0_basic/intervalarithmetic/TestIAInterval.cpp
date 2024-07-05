#include <vector>
#include "TestEnvironment.hpp"

using namespace ::testing;
using std::string;
using std::vector;

struct IAIntervalTestParameter {
  double infA;
  double supA;
  double infB;
  double supB;
  double infAB;
  double supAB;
};

class IAIntervalTest : public TestWithParam<IAIntervalTestParameter> {
protected:
  IAInterval a, b, ab;
public:
  IAIntervalTest() {
    a = IAInterval(GetParam().infA, GetParam().supA);
    b = IAInterval(GetParam().infB, GetParam().supB);
    ab = IAInterval(GetParam().infAB, GetParam().supAB);
  }
};

// Summation tests
typedef IAIntervalTest SummationIAIntervalTest;

TEST_P(SummationIAIntervalTest, AdditionTest) { EXPECT_EQ(a + b, ab); }

INSTANTIATE_TEST_CASE_P(IAIntervalTest, SummationIAIntervalTest,
                        Values(IAIntervalTestParameter{1.0, 2.0, 3.0, 4.0, 4.0, 6.0},
                               IAIntervalTestParameter{1.0, 2.0, -3.0, 4.0, -2.0, 6.0},
                               IAIntervalTestParameter{1.0, 2.0, -4.0, -3.0, -3.0, -1.0},
                               IAIntervalTestParameter{-1.0, 2.0, 3.0, 4.0, 2.0, 6.0},
                               IAIntervalTestParameter{-1.0, 2.0, -3.0, 4.0, -4.0, 6.0},
                               IAIntervalTestParameter{-1.0, 2.0, -4.0, -3.0, -5.0, -1.0},
                               IAIntervalTestParameter{-2.0, -1.0, 3.0, 4.0, 1.0, 3.0},
                               IAIntervalTestParameter{-2.0, -1.0, -3.0, 4.0, -5.0, 3.0},
                               IAIntervalTestParameter{-2.0, -1.0, -4.0, -3.0, -6.0, -4.0}));

// Difference tests
typedef IAIntervalTest DifferenceIAIntervalTest;

TEST_P(DifferenceIAIntervalTest, DifferenceTest) { EXPECT_EQ(a - b, ab); }

INSTANTIATE_TEST_CASE_P(IAIntervalTest, DifferenceIAIntervalTest,
                        Values(IAIntervalTestParameter{1.0, 2.0, 3.0, 4.0, -3.0, -1.0},
                               IAIntervalTestParameter{1.0, 2.0, -3.0, 4.0, -3.0, 5.0},
                               IAIntervalTestParameter{1.0, 2.0, -4.0, -3.0, 4.0, 6.0},
                               IAIntervalTestParameter{-1.0, 2.0, 3.0, 4.0, -5.0, -1.0},
                               IAIntervalTestParameter{-1.0, 2.0, -3.0, 4.0, -5.0, 5.0},
                               IAIntervalTestParameter{-1.0, 2.0, -4.0, -3.0, 2.0, 6.0},
                               IAIntervalTestParameter{-2.0, -1.0, 3.0, 4.0, -6.0, -4.0},
                               IAIntervalTestParameter{-2.0, -1.0, -3.0, 4.0, -6.0, 2.0},
                               IAIntervalTestParameter{-2.0, -1.0, -4.0, -3.0, 1.0, 3.0}));

// Multiplication tests
typedef IAIntervalTest MultiplicationIAIntervalTest;

TEST_P(MultiplicationIAIntervalTest, MultiplicationTest) { EXPECT_EQ(a * b, ab); }

INSTANTIATE_TEST_CASE_P(IAIntervalTest, MultiplicationIAIntervalTest,
                        Values(IAIntervalTestParameter{1.0, 2.0, 3.0, 4.0, 3.0, 8.0},
                               IAIntervalTestParameter{1.0, 2.0, -3.0, 4.0, -6.0, 8.0},
                               IAIntervalTestParameter{1.0, 2.0, -4.0, -3.0, -8.0, -3.0},
                               IAIntervalTestParameter{-1.0, 2.0, 3.0, 4.0, -4.0, 8.0},
                               IAIntervalTestParameter{-1.0, 2.0, -3.0, 4.0, -6.0, 8.0},
                               IAIntervalTestParameter{-1.0, 2.0, -4.0, -3.0, -8.0, 4.0},
                               IAIntervalTestParameter{-2.0, -1.0, 3.0, 4.0, -8.0, -3.0},
                               IAIntervalTestParameter{-2.0, -1.0, -3.0, 4.0, -8.0, 6.0},
                               IAIntervalTestParameter{-2.0, -1.0, -4.0, -3.0, 3.0, 8.0}));

// Division tests
typedef IAIntervalTest DivisionIAIntervalTest;

TEST_P(DivisionIAIntervalTest, DivisionTest) { EXPECT_EQ(a / b, ab); }

INSTANTIATE_TEST_CASE_P(IAIntervalTest, DivisionIAIntervalTest,
                        Values(IAIntervalTestParameter{4.0, 8.0, 1.0, 2.0, 2.0, 8.0},
                               IAIntervalTestParameter{4.0, 8.0, -2.0, -1.0, -8.0, -2.0},
                               IAIntervalTestParameter{-4.0, 8.0, 1.0, 2.0, -4.0, 8.0},
                               IAIntervalTestParameter{-4.0, 8.0, -2.0, -1.0, -8.0, 4.0},
                               IAIntervalTestParameter{-8.0, -4.0, 1.0, 2.0, -8.0, -2.0},
                               IAIntervalTestParameter{-8.0, -4.0, -2.0, -1.0, 2.0, 8.0}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
