#include <complex>
#include <utility>
#include "SaveLoad.hpp"
#include "TestEnvironment.hpp"

template<typename T>
T testSaveLoad(T input) {
  Saver saver("SaveLoadTest");
  saver << input;
  saver.close();
  Loader loader("SaveLoadTest");
  T output;
  loader >> output;
  loader.close();
  return output;
}

TEST(SaveLoadTest, IAIntervalTest) {
  double lb = double(rand()) / RAND_MAX;
  double ub = lb + double(rand()) / RAND_MAX;
  IAInterval input(lb, ub);
  IAInterval output = testSaveLoad(input);
  ASSERT_EQ(input, output);
}

TEST(SaveLoadTest, IACIntervalTest) {
  double lb = double(rand()) / RAND_MAX;
  double ub = lb + double(rand()) / RAND_MAX;
  double im_lb = double(rand()) / RAND_MAX;
  double im_ub = im_lb + double(rand()) / RAND_MAX;
  IACInterval input(IAInterval(lb, ub), IAInterval(im_lb, im_ub));
  IACInterval output = testSaveLoad(input);
  ASSERT_EQ(input, output);
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
