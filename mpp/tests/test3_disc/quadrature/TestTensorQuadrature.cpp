#include <vector>
#include "Quadrature.hpp"
#include "QuadratureGL.hpp"
#include "TestEnvironment.hpp"

using namespace ::testing;
using std::initializer_list;
using std::vector;

struct TensorQuadratureTestParameter {
  std::vector<QuadratureT<double, SpaceDimension, TimeDimension>> quads;
  vector<double> expectedWeights;
  vector<Point> exptectedPoints;
};

class TensorQuadratureTest : public TestWithParam<TensorQuadratureTestParameter> {
public:
  QTensor quad;
  vector<double> expectedWeights;
  vector<Point> expectedPoints;
  int tensor_dim;

  TensorQuadratureTest() :
      quad(GetParam().quads), expectedWeights(GetParam().expectedWeights),
      expectedPoints(GetParam().exptectedPoints), tensor_dim(GetParam().quads.size()){


                                                  };
};

TEST_P(TensorQuadratureTest, QuadratureTest) {
  for (int n = 0; n < quad.size(); ++n) {
    EXPECT_NEAR(expectedWeights[n], quad.Weight(n), 0);
    for (int d = 0; d < tensor_dim; d++) {
      EXPECT_NEAR(expectedPoints[n][d], quad.QPoint(n)[d], 0);
    }
  }
}

INSTANTIATE_TEST_SUITE_P(TwoDimTensorQuadratureTest, TensorQuadratureTest,
                         Values(TensorQuadratureTestParameter{
                             std::vector<Quadrature>{Qint(2), Qint(2)},
                             {0.25, 0.25, 0.25, 0.25},
                             {Point(0.2113248654051871177454256097490212721761991243649365619906,
                                    0.2113248654051871177454256097490212721761991243649365619906),
                              Point(0.7886751345948128822545743902509787278238008756350634380093,
                                    0.2113248654051871177454256097490212721761991243649365619906),
                              Point(0.2113248654051871177454256097490212721761991243649365619906,
                                    0.7886751345948128822545743902509787278238008756350634380093),
                              Point(0.7886751345948128822545743902509787278238008756350634380093,
                                    0.7886751345948128822545743902509787278238008756350634380093)}})

);

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}