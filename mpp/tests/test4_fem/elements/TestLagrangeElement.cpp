#include "LagrangeDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "ScalarElement.hpp"
#include "TestEnvironment.hpp"
#include "Vector.hpp"

using std::vector;

const double testTol = 1e-12;

struct TestParams {
  std::string meshesName;
  int degree;
};

static vector<TestParams> testCasesInt{TestParams{"Interval", 0}, TestParams{"Interval", 1},
                                       TestParams{"Interval", 2}, TestParams{"Interval", 3},
                                       TestParams{"Interval", 4}, TestParams{"Interval", 5},
                                       TestParams{"Interval", 6}, TestParams{"Interval", 7},
                                       TestParams{"Interval", 8}};

static vector<TestParams>
    testCasesTri{TestParams{"Triangle", 0},         TestParams{"Triangle", 1},
                 TestParams{"Triangle", 2},         TestParams{"Triangle", 3},
                 TestParams{"Triangle", 4},         TestParams{"Triangle", 5},
                 TestParams{"Triangle", 6},         TestParams{"Square2Triangles", 0},
                 TestParams{"Square2Triangles", 1}, TestParams{"Square2Triangles", 2},
                 TestParams{"Square2Triangles", 3}, TestParams{"Square2Triangles", 4},
                 TestParams{"Square2Triangles", 5}, TestParams{"Square2Triangles", 6}};

static vector<TestParams> testCasesQuad{TestParams{"Square", 0}, TestParams{"Square", 1},
                                        TestParams{"Square", 2}, TestParams{"Square", 3},
                                        TestParams{"Square", 4}, TestParams{"Square", 5},
                                        TestParams{"Square", 6}, TestParams{"Square", 7},
                                        TestParams{"Square", 8}};

static vector<TestParams> testCasesTet{
    TestParams{"Tetrahedron", 0}, TestParams{"Tetrahedron", 1}, TestParams{"Tetrahedron", 2},
    TestParams{"Tetrahedron", 3} //,
    //    TestParams{"Tetrahedron", 4}
};

static vector<TestParams> testCasesHex{TestParams{"Hexahedron", 0}, TestParams{"Hexahedron", 1},
                                       TestParams{"Hexahedron", 2}, TestParams{"Hexahedron", 3},
                                       TestParams{"Hexahedron", 4}, TestParams{"Hexahedron", 5},
                                       TestParams{"Hexahedron", 6}, TestParams{"Hexahedron", 7},
                                       TestParams{"Hexahedron", 8}};

/// To avoid the same code in each test
#define LAGRANGE_TESTS(testclass, testCases)                                                       \
                                                                                                   \
  TEST_P(testclass, TestContinuity) { checkContinuity(); }                                         \
                                                                                                   \
  INSTANTIATE_TEST_CASE_P(TestLagrangeElement, testclass, ValuesIn(testCases));

class TestLagrangeElement : public TestWithParam<TestParams> {
protected:
  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const LagrangeDiscretization> disc;
  Vector *vec = nullptr;

  TestLagrangeElement() {
    meshes = MeshesCreator(GetParam().meshesName).WithPLevel(2).WithLevel(2).CreateShared();
    disc = std::make_shared<const LagrangeDiscretization>(*meshes, GetParam().degree);
    vec = new Vector(0.0, disc);
  }

  ~TestLagrangeElement() { delete vec; }

  virtual void checkContinuity() = 0;
};

class IntervalMeshTest : public TestLagrangeElement {
protected:
  void checkContinuity() override {
    if (GetParam().degree == 0) return;
    const Mesh &mesh = meshes->fine();
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      ScalarElement E_c(*vec, *c);
      for (int face = 0; face < c.Faces(); ++face) {
        if (E_c.Bnd(face) != -1) continue;
        Point p = c.FaceCorner(face, 0);
        Point p_c = c.GlobalToLocal(p);
        cell neighbour = vec->find_neighbour_cell(c, face);
        Point p_n = neighbour.GlobalToLocal(p);
        ScalarElement E_n(*vec, *neighbour);

        Point np = c.FaceCorner(face, 0);
        int i_c = E_c.NodalPointIndex(np);
        int i_n = E_n.NodalPointIndex(np);
        EXPECT_NEAR(E_c.Value(p_c, i_c), E_n.Value(p_n, i_n), testTol);
      }
    }
  }
};

LAGRANGE_TESTS(IntervalMeshTest, testCasesInt)

#if SpaceDimension >= 2

class TestLagrangeElement2D : public TestLagrangeElement {
  Point GlobalPointOnFace(const cell &c, int face, double t) {
    return (*c).FaceLocalToGlobal(face, Point(t, 0.0));
  }
protected:
  void checkContinuity() override {
    if (GetParam().degree == 0) return;
    for (int cnt = 0; cnt < 10; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        ScalarElement E_c(*vec, *c);
        for (int face = 0; face < c.Faces(); ++face) {
          if (E_c.Bnd(face) != -1) continue;
          Point p = GlobalPointOnFace(c, face, RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face);
          Point p_n = neighbour.GlobalToLocal(p);
          ScalarElement E_n(*vec, *neighbour);

          for (int i = 0; i <= GetParam().degree; ++i) {
            Point np = GlobalPointOnFace(c, face, double(i) / GetParam().degree);
            int i_c = E_c.NodalPointIndex(np);
            int i_n = E_n.NodalPointIndex(np);
            EXPECT_NEAR(E_c.Value(p_c, i_c), E_n.Value(p_n, i_n), testTol);
          }
        }
      }
    }
  }
};

using TriangularMeshTest = TestLagrangeElement2D;

LAGRANGE_TESTS(TriangularMeshTest, testCasesTri)

using QuadraticMeshTest = TestLagrangeElement2D;

LAGRANGE_TESTS(QuadraticMeshTest, testCasesQuad)

#endif
#if SpaceDimension >= 3

class TetrahedronMeshTest : public TestLagrangeElement {
  Point GlobalPointOnFace(const cell &c, int face, double t, double s) {
    return (*c).FaceLocalToGlobal(face, Point(t, s));
  }
protected:
  void checkContinuity() override {
    if (GetParam().degree == 0) return;
    for (int cnt = 0; cnt < 10; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        ScalarElement E_c(*vec, *c);
        for (int face = 0; face < c.Faces(); ++face) {
          if (E_c.Bnd(face) != -1) continue;
          double t = RandomDouble(0.0, 1.0);
          Point p = GlobalPointOnFace(c, face, t, RandomDouble(0.0, t));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face);
          Point p_n = neighbour.GlobalToLocal(p);
          ScalarElement E_n(*vec, *neighbour);

          for (int i = 0; i <= GetParam().degree; ++i) {
            for (int j = 0; j <= GetParam().degree - i; ++j) {
              Point np = GlobalPointOnFace(c, face, double(i) / GetParam().degree,
                                           double(j) / GetParam().degree);
              int i_c = E_c.NodalPointIndex(np);
              int i_n = E_n.NodalPointIndex(np);
              EXPECT_NEAR(E_c.Value(p_c, i_c), E_n.Value(p_n, i_n), testTol);
            }
          }
        }
      }
    }
  }
};

LAGRANGE_TESTS(TetrahedronMeshTest, testCasesTet)

class HexahedronMeshTest : public TestLagrangeElement {
  Point GlobalPointOnFace(const cell &c, int face, double t, double s) {
    return (*c).FaceLocalToGlobal(face, Point(t, s));
  }
protected:
  void checkContinuity() override {
    if (GetParam().degree == 0) return;
    for (int cnt = 0; cnt < 5; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        ScalarElement E_c(*vec, *c);
        for (int face = 0; face < c.Faces(); ++face) {
          if (E_c.Bnd(face) != -1) continue;
          Point p = GlobalPointOnFace(c, face, RandomDouble(0.0, 1.0), RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face);
          Point p_n = neighbour.GlobalToLocal(p);
          ScalarElement E_n(*vec, *neighbour);

          for (int i = 0; i <= GetParam().degree; ++i) {
            for (int j = 0; j <= GetParam().degree; ++j) {
              Point np = GlobalPointOnFace(c, face, double(i) / GetParam().degree,
                                           double(j) / GetParam().degree);
              int i_c = E_c.NodalPointIndex(np);
              int i_n = E_n.NodalPointIndex(np);
              EXPECT_NEAR(E_c.Value(p_c, i_c), E_n.Value(p_n, i_n), testTol);
            }
          }
        }
      }
    }
  }
};

LAGRANGE_TESTS(HexahedronMeshTest, testCasesHex)

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
