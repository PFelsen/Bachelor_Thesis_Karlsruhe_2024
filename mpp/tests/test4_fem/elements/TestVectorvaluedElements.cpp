#include "LagrangeDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "TensorElement.hpp"
#include "TestEnvironment.hpp"
#include "Vector.hpp"
#include "VectorFieldElement.hpp"

#include <vector>

using std::vector;

const double ansatzTol = 1e-12;

class TestVectorvaluedElement : public TestWithParam<std::pair<std::string, int>> {
protected:
  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const LagrangeDiscretization> disc;
  std::shared_ptr<Vector> vec;

  TestVectorvaluedElement() {
    meshes = MeshesCreator(GetParam().first).WithPLevel(2).WithLevel(2).CreateShared();
    disc = std::make_shared<LagrangeDiscretization>(*meshes, 1, GetParam().second);
    vec = std::make_shared<Vector>(0.0, disc);
  }
};

TEST_P(TestVectorvaluedElement, TestVectorFieldValue) {
  VectorFieldElement element(*vec, *vec->cells());
  for (int q = 0; q < element.nQ(); ++q) {
    for (int i = 0; i < element.size(); ++i) {
      for (int k = 0; k < GetParam().second; ++k) {
        VectorFieldComponent v = element.VectorComponentValue(q, i, k);
        VectorField V(k, element.Value(q, i));
        EXPECT_EQ(VectorField(v), V);
      }
    }
  }
}

TEST_P(TestVectorvaluedElement, TestVectorFieldDerivative) {
  VectorFieldElement element(*vec, *vec->cells());
  for (int q = 0; q < element.nQ(); ++q) {
    for (int i = 0; i < element.size(); ++i) {
      for (int k = 0; k < GetParam().second; ++k) {
        TensorRow t = element.VectorRowGradient(q, i, k);
        Tensor T(k, element.Derivative(q, i));
        EXPECT_EQ(Tensor(t), T);
      }
    }
  }
}

TEST_P(TestVectorvaluedElement, TestTensorValue) {
  TensorElement element(*vec, *vec->cells());
  for (int q = 0; q < element.nQ(); ++q) {

    for (int i = 0; i < element.size(); ++i) {
      for (int k = 0; k < GetParam().second; ++k) {
        TensorComponent c = element.TensorComponentValue(q, i, k);
        Tensor C;
        if (k < element.Dim()) C = Tensor(0, element.VectorComponentValue(q, i, k));
        else if (k >= 2 * element.Dim())
          C = Tensor(2, element.VectorComponentValue(q, i, k - 2 * element.Dim()));
        else C = Tensor(1, element.VectorComponentValue(q, i, k - element.Dim()));
        EXPECT_EQ(Tensor(c), C);
      }
    }
  }
}

static std::vector<std::pair<std::string, int>> testCases{
#if SpaceDimension >= 2
    std::make_pair("Triangle", 1),
    std::make_pair("Triangle", 2),
    std::make_pair("Square2Triangles", 1),
    std::make_pair("Square2Triangles", 2),
    std::make_pair("Square", 1),
    std::make_pair("Square", 2),
#endif
#if SpaceDimension >= 3
    std::make_pair("Tetrahedron", 1),
    std::make_pair("Tetrahedron", 2),
    std::make_pair("Hexahedron", 1),
    std::make_pair("Hexahedron", 2)
#endif
};

INSTANTIATE_TEST_SUITE_P(TestVectorvaluedElement, TestVectorvaluedElement, ValuesIn(testCases));

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(TestVectorvaluedElement);

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}