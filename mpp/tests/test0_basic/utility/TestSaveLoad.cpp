#include <complex>
#include <utility>
#include "TestEnvironment.hpp"
#include "utility/SaveLoad.hpp"

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

TEST(SaveLoadTest, ShortIntegerTest) {
  short int input = rand();
  short int output = testSaveLoad(input);
  ASSERT_EQ(input, output);
}

TEST(SaveLoadTest, IntegerTest) {
  int input = rand();
  int output = testSaveLoad(input);
  ASSERT_EQ(input, output);
}

TEST(SaveLoadTest, Size_tTest) {
  size_t input = rand();
  size_t output = testSaveLoad(input);
  ASSERT_EQ(input, output);
}

TEST(SaveLoadTest, DoubleTest) {
  double input = double(rand()) / RAND_MAX;
  double output = testSaveLoad(input);
  ASSERT_DOUBLE_EQ(input, output);
}

TEST(SaveLoadTest, CharTest) {
  char input = rand();
  char output = testSaveLoad(input);
  ASSERT_EQ(input, output);
}

TEST(SaveLoadTest, ComplexTest) {
  double re = double(rand()) / RAND_MAX;
  double im = double(rand()) / RAND_MAX;
  std::complex<double> input(re, im);
  std::complex<double> output = testSaveLoad(input);
  ASSERT_DOUBLE_EQ(re, std::real(output));
  ASSERT_DOUBLE_EQ(im, std::imag(output));
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
