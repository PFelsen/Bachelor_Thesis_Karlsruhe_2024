#include "TestEnvironment.hpp"

void function(int p) {
  if (p == 0 || p == 2) THROW("test error")
}

void testTryCatch() {
  TRY {
    std::vector<double> v(1);
    if (PPM->Master(0)) v.at(2);
  }
  CATCH("Error on proc" + std::to_string(PPM->Proc(0)))
}

void testThrow(){TRY{function(PPM -> Proc(0));
}
CATCH("")
}

/// Tests do not work in parallel yet
TEST(AssertionTest, ParallelTryCatchTest) {
  EXPECT_EXIT(testTryCatch(), testing::ExitedWithCode(1), "");
}

TEST(AssertionTest, ParallelThrowTest) { EXPECT_EXIT(testThrow(), testing::ExitedWithCode(1), ""); }

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
