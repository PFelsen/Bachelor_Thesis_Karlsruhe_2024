#include "TestTensorInvert.hpp"
#include "Tensor.hpp"
#include "TestEnvironment.hpp"

template<int sDim>
class TensorInvertTest : public TestWithParam<TensorData> {
protected:
  using TENSOR = TensorT<double, sDim>;

  TENSOR T;
  TENSOR invT;
  TENSOR transposeInversT;

  TensorInvertTest() {
    mpp_geometry::SetTolerance(1e-14);
    for (int i = 0; i < sDim; ++i)
      for (int j = 0; j < sDim; ++j) {
        T[i][j] = GetParam().tensorValues[3 * i + j];
        invT[i][j] = GetParam().inverseValues[3 * i + j];
        transposeInversT[i][j] = GetParam().inverseValues[3 * j + i];
      }
  }
};

/// To avoid the same code in each test
#define TENSOR_INVERT_TESTS(tensorClass)                                                           \
                                                                                                   \
  TEST_P(tensorClass, InvertTest) { EXPECT_EQ(Invert(T), invT); }                                  \
                                                                                                   \
  TEST_P(tensorClass, TransposeInvertTest) { EXPECT_EQ(transposeInvert(T), transposeInversT); }

using TensorInvertTest_1 = TensorInvertTest<1>;

TENSOR_INVERT_TESTS(TensorInvertTest_1)

INSTANTIATE_TEST_CASE_P(TensorInvertTest, TensorInvertTest_1, ValuesIn(tensorInvertData_1));

using TensorInvertTest_2 = TensorInvertTest<2>;

TENSOR_INVERT_TESTS(TensorInvertTest_2)

INSTANTIATE_TEST_CASE_P(TensorInvertTest, TensorInvertTest_2,
                        ValuesIn(combineTensorData(tensorInvertData_1, tensorInvertData_2)));

using TensorInvertTest_3 = TensorInvertTest<3>;

TENSOR_INVERT_TESTS(TensorInvertTest_3)

INSTANTIATE_TEST_CASE_P(TensorInvertTest, TensorInvertTest_3,
                        ValuesIn(combineTensorData(tensorInvertData_1, tensorInvertData_2,
                                                   tensorInvertData_3)));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
