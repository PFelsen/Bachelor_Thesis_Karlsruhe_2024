#include "LagrangeNodalPoints.hpp"
#include "MeshesCreator.hpp"
#include "RTDiscretization.hpp"
#include "RTElement.hpp"
#include "TestEnvironment.hpp"
#include "Vector.hpp"

#include <vector>

using std::vector;

const double testTol = 5e-12;

struct TestParams {
  std::string meshesName;
  int order;
};

static vector<TestParams> testCasesInt{TestParams{"Interval", 0}};

static vector<TestParams>
    testCasesTri{TestParams{"Triangle", 0},         TestParams{"Triangle", 1},
                 TestParams{"Triangle", 2},         TestParams{"Triangle", 3},
                 TestParams{"Triangle", 4},         TestParams{"Triangle", 5},
                 TestParams{"Triangle", 6},         TestParams{"Triangle", 7},
                 TestParams{"Square2Triangles", 0}, TestParams{"Square2Triangles", 1},
                 TestParams{"Square2Triangles", 2}, TestParams{"Square2Triangles", 3},
                 TestParams{"Square2Triangles", 4}, TestParams{"Square2Triangles", 5},
                 TestParams{"Square2Triangles", 6}, TestParams{"Square2Triangles", 7}};

static vector<TestParams> testCasesQuad{TestParams{"Square", 0}};

static vector<TestParams> testCasesHex{TestParams{"Hexahedron", 0}};

/// To avoid the same code in each test
#define RT_TESTS(testclass, testCases)                                                             \
                                                                                                   \
  TEST_P(testclass, TestContinuity) { checkContinuity(); }                                         \
                                                                                                   \
  TEST_P(testclass, TestNormalConsistency) { checkNormal(); }                                      \
                                                                                                   \
  TEST_P(testclass, TestMaxK) { checkMaxK(); }                                                     \
                                                                                                   \
  TEST_P(testclass, TestIndexingK) { checkIndexingK(); }                                           \
                                                                                                   \
  INSTANTIATE_TEST_CASE_P(TestRTElement, testclass, ValuesIn(testCases));

class TestRTElement : public TestWithParam<TestParams> {
protected:
  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const RTDiscretization> disc;
  std::shared_ptr<Vector> vec;

  TestRTElement() {
    meshes = MeshesCreator(GetParam().meshesName).WithPLevel(2).WithLevel(2).CreateShared();
    disc = std::make_shared<const RTDiscretization>(*meshes, GetParam().order);
    vec = std::make_shared<Vector>(0.0, disc);
  }

  virtual void checkContinuity() = 0;

  void checkNormal() {
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      RTElement E_c(*vec, *c);
      for (int face_c = 0; face_c < c.Faces(); ++face_c) {
        if (E_c.Bnd(face_c) != -1) continue;
        cell neighbour = vec->find_neighbour_cell(c, face_c);
        RTElement E_n(*vec, *neighbour);
        int face_n = (*neighbour).Face(c.Face(face_c));

        EXPECT_EQ(E_c.OrientedNormal(face_c), E_n.OrientedNormal(face_n));
      }
    }
  }

  void checkMaxK() {
    RTElement element(*vec, *vec->cells());
    int i = 0;
    for (; i < vec->cells().Faces() * (GetParam().order + 1); ++i)
      EXPECT_EQ(element.get_maxk(i), 1);
    for (; i < element.size(); ++i)
      EXPECT_EQ(element.get_maxk(i), 2);
  }

  void checkIndexingK() {
    RTElement element(*vec, *vec->cells());
    for (int m = 0; m < 3; ++m)
      for (int i = 0; i < element.size(); ++i)
        for (int k = 0; k < element.get_maxk(i); ++k)
          EXPECT_EQ(element.indexing_k(i, k, m), element.get_maxk(i) * m + k);
  }
};

class IntervalMeshTest : public TestRTElement {
protected:
  void checkContinuity() override {
    const Mesh &mesh = meshes->fine();
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      RTElement E_c(*vec, *c);
      for (int face = 0; face < c.Faces(); ++face) {
        if (E_c.Bnd(face) != -1) continue;
        Point p = c.FaceCorner(face, 0);
        Point p_c = c.GlobalToLocal(p);
        cell neighbour = vec->find_neighbour_cell(c, face);
        Point p_n = neighbour.GlobalToLocal(p);
        RTElement E_n(*vec, *neighbour);

        Point np = c.FaceCorner(face, 0);
        int i_c = E_c.NodalPointIndex(np);
        int i_n = E_n.NodalPointIndex(np);
        EXPECT_NEAR(E_c.VelocityField(p_c, i_c) * E_c.OrientedNormal(face),
                    E_n.VelocityField(p_n, i_n) * E_n.OrientedNormal(face), testTol);
      }
    }
  }
};

RT_TESTS(IntervalMeshTest, testCasesInt)

#if SpaceDimension >= 2

class TestRTElement2D : public TestRTElement {
  Point GlobalPointOnFace(const cell &c, int face, double t) {
    return (*c).FaceLocalToGlobal(face, Point(t, 0.0));
  }
protected:
  void checkContinuity() override {
    for (int cnt = 0; cnt < 10; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        RTElement E_c(*vec, *c);
        for (int face_c = 0; face_c < c.Faces(); ++face_c) {
          if (E_c.Bnd(face_c) != -1) continue;
          Point p = GlobalPointOnFace(c, face_c, RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face_c);
          Point p_n = neighbour.GlobalToLocal(p);
          RTElement E_n(*vec, *neighbour);
          int face_n = (*neighbour).Face(c.Face(face_c));

          for (int i = 1; i <= GetParam().order + 1; ++i) {
            Point np = GlobalPointOnFace(c, face_c, double(i) / (GetParam().order + 2));
            int i_c = E_c.NodalPointIndex(np);
            int i_n = E_n.NodalPointIndex(np);
            EXPECT_NEAR(E_c.VelocityField(p_c, i_c) * E_c.OrientedNormal(face_c),
                        E_n.VelocityField(p_n, i_n) * E_n.OrientedNormal(face_n), testTol);
          }
        }
      }
    }
  }
};

using TriangularMeshTest = TestRTElement2D;

RT_TESTS(TriangularMeshTest, testCasesTri)

using QuadraticMeshTest = TestRTElement2D;

RT_TESTS(QuadraticMeshTest, testCasesQuad)

#endif
#if SpaceDimension >= 3

class HexahedronMeshTest : public TestRTElement {
  Point GlobalPointOnFace(const cell &c, int face, double t, double s) {
    return (*c).FaceLocalToGlobal(face, Point(t, s));
  }
protected:
  void checkContinuity() override {
    for (int cnt = 0; cnt < 10; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        RTElement E_c(*vec, *c);
        for (int face_c = 0; face_c < c.Faces(); ++face_c) {
          if (E_c.Bnd(face_c) != -1) continue;
          Point p = GlobalPointOnFace(c, face_c, RandomDouble(0.0, 1.0), RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face_c);
          Point p_n = neighbour.GlobalToLocal(p);
          RTElement E_n(*vec, *neighbour);
          int face_n = (*neighbour).Face(c.Face(face_c));

          Point np = c.Face(face_c);
          int i_c = E_c.NodalPointIndex(np);
          int i_n = E_n.NodalPointIndex(np);
          EXPECT_NEAR(E_c.VelocityField(p_c, i_c) * E_c.OrientedNormal(face_c),
                      E_n.VelocityField(p_n, i_n) * E_n.OrientedNormal(face_n), testTol);
        }
      }
    }
  }
};

RT_TESTS(HexahedronMeshTest, testCasesHex)

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
