#include "QuadratureGL.hpp"

#include "TestEnvironment.hpp"

struct GLQuadTestParam {
  Quadrature quad;
  int exactUpTo;

  friend std::ostream &operator<<(std::ostream &os, const GLQuadTestParam &testparam) {
    return os << testparam.quad.Name() << " exact up to degree: " << testparam.exactUpTo;
  }
};

class TestGaussLobattoQuadrature : public TestWithParam<GLQuadTestParam> {
protected:
  Quadrature quad;
  std::vector<double> calculated_values;

  TestGaussLobattoQuadrature() :
      quad(GetParam().quad), calculated_values(GetParam().exactUpTo + 2) {
    for (int d = 0; d < calculated_values.size(); ++d) {
      calculated_values[d] = evalQuadrature(d);
    }
  }

  double evalQuadrature(int exponent) {
    double val = 0.0;
    for (int q = 0; q < quad.size(); ++q)
      val += quad.Weight(q) * pow(quad.QPoint(q)[0], exponent);
    return val;
  }
};

TEST_P(TestGaussLobattoQuadrature, TestGaussLobattoQuadrature) {
  for (int d = 0; d <= GetParam().exactUpTo; ++d) {
    EXPECT_NEAR(1.0 / (d + 1), calculated_values[d], 1e-15);
  }
  int fail_degree = GetParam().exactUpTo + 1;

  EXPECT_GT(abs(1.0 / (fail_degree + 1) - calculated_values[fail_degree]), 1e-15);
}

INSTANTIATE_TEST_SUITE_P(QuadratureTest, TestGaussLobattoQuadrature,
                         Values(GLQuadTestParam{QintGaussLobatto(1), 1},
                                GLQuadTestParam{QintGaussLobatto(2), 2 * 2 - 3},
                                GLQuadTestParam{QintGaussLobatto(3), 2 * 3 - 3},
                                GLQuadTestParam{QintGaussLobatto(4), 2 * 4 - 3},
                                GLQuadTestParam{QintGaussLobatto(5), 2 * 5 - 3},
                                GLQuadTestParam{QintGaussLobatto(6), 2 * 6 - 3},
                                GLQuadTestParam{QintGaussLobatto(7), 2 * 7 - 3},
                                GLQuadTestParam{QintGaussLobatto(8), 2 * 8 - 3},
                                GLQuadTestParam{QintGaussLobatto(9), 2 * 9 - 3},
                                GLQuadTestParam{QintGaussLobatto(10), 2 * 10 - 3},
                                GLQuadTestParam{QintGaussLobatto(11), 2 * 11 - 3},
                                GLQuadTestParam{QintGaussLobatto(12), 2 * 12 - 3},
                                GLQuadTestParam{QintGaussLobatto(13), 2 * 13 - 3}));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}