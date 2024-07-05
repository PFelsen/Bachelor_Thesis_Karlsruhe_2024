#include "ArgyrisDiscretization.hpp"
#include "ArgyrisElement.hpp"
#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"
#include "Vector.hpp"

#include <vector>

using std::vector;

const double testTol = 1e-12;

#if SpaceDimension >= 2

class TestArgyrisElement : public TestWithParam<std::string> {
  Point GlobalPointOnFace(const cell &c, int face, double t) {
    return (*c).FaceLocalToGlobal(face, Point(t, 0.0));
  }
protected:
  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const ArgyrisDiscretization> disc;
  std::shared_ptr<Vector> vec;

  TestArgyrisElement() {
    meshes = MeshesCreator(GetParam()).WithPLevel(2).WithLevel(2).CreateShared();
    disc = std::make_shared<const ArgyrisDiscretization>(*meshes);
    vec = std::make_shared<Vector>(0.0, disc);
  }

  // Placed here to exploit friendship
  void checkIndexingShape() {
    ArgyrisElement element(*vec, *vec->cells());
    for (int i = 0; i < element.size(); ++i)
      for (int k = 0; k < element.get_maxk(i); ++k) {
        if (i < 3) {
          EXPECT_EQ(element.indexingShape(i, k), i * 6 + k);
        } else {
          EXPECT_EQ(element.indexingShape(i, k), 15 + i);
        }
      }
  }

  // Placed here to exploit friendship
  void checkNodalPointValues() {
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      ArgyrisElement element(*vec, *c);

      for (int corner = 0; corner < c.Corners(); ++corner) {
        for (int i = 0; i < element.size(); ++i)
          for (int k = 0; k < element.get_maxk(i); ++k) {
            int j = element.indexingShape(i, k);
            if (j == corner * 6) {
              EXPECT_NEAR(element.Value(c.LocalCorner(corner), i, k), 1.0, testTol);
            } else {
              EXPECT_NEAR(element.Value(c.LocalCorner(corner), i, k), 0.0, testTol);
            }
            VectorField d = element.Derivative(c.LocalCorner(corner), i, k);
            for (int l = 0; l < 2; ++l) {
              if (j == corner * 6 + l + 1) {
                EXPECT_NEAR(d[l], 1.0, testTol);
              } else {
                EXPECT_NEAR(d[l], 0.0, testTol);
              }
            }
            SymTensor h = element.Hessian(c.LocalCorner(corner), i, k);
            if (j == corner * 6 + 3) {
              EXPECT_NEAR(h(0, 0), 1.0, testTol);
            } else {
              EXPECT_NEAR(h(0, 0), 0.0, testTol);
            }
            if (j == corner * 6 + 4) {
              EXPECT_NEAR(h(0, 1), 1.0, testTol);
            } else {
              EXPECT_NEAR(h(0, 1), 0.0, testTol);
            }
            if (j == corner * 6 + 5) {
              EXPECT_NEAR(h(1, 1), 1.0, testTol);
            } else {
              EXPECT_NEAR(h(1, 1), 0.0, testTol);
            }
          }
      }

      for (int face = 0; face < c.Faces(); ++face) {
        for (int i = 0; i < element.size(); ++i)
          for (int k = 0; k < element.get_maxk(i); ++k) {
            int j = element.indexingShape(i, k);
            double value =
                element.Derivative(c.LocalFace(face), i, k) * element.OrientedNormal(face);
            if (j == 18 + face) {
              EXPECT_NEAR(value, 1.0, testTol);
            } else {
              EXPECT_NEAR(value, 0.0, testTol);
            }
          }
      }
    }
  }

  void checkContinuity() {
    for (int cnt = 0; cnt < 10; ++cnt) {
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        ArgyrisElement E_c(*vec, *c);
        for (int face = 0; face < c.Faces(); ++face) {
          if (E_c.Bnd(face) != -1) continue;
          Point p = GlobalPointOnFace(c, face, RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face);
          Point p_n = neighbour.GlobalToLocal(p);
          ArgyrisElement E_n(*vec, *neighbour);

          // shape functions corresponding to left face corner
          Point np = GlobalPointOnFace(c, face, 0.0);
          int i_c = E_c.NodalPointIndex(np);
          int i_n = E_n.NodalPointIndex(np);
          for (int k = 0; k < 6; ++k)
            EXPECT_NEAR(E_c.Value(p_c, 6 * i_c + k), E_n.Value(p_n, 6 * i_n + k), testTol);

          // shape functions corresponding to right face corner
          np = GlobalPointOnFace(c, face, 1.0);
          i_c = E_c.NodalPointIndex(np);
          i_n = E_n.NodalPointIndex(np);
          for (int k = 0; k < 6; ++k)
            EXPECT_NEAR(E_c.Value(p_c, 6 * i_c + k), E_n.Value(p_n, 6 * i_n + k), testTol);

          // shape functions corresponding to face middle
          np = GlobalPointOnFace(c, face, 0.5);
          i_c = E_c.NodalPointIndex(np);
          i_n = E_n.NodalPointIndex(np);
          EXPECT_NEAR(E_c.Value(p_c, 18 + i_c), E_n.Value(p_n, 18 + i_n), testTol);
        }
      }
    }
  }

  void checkDerivativeContinuity() {
    for (int cnt = 0; cnt < 10; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        ArgyrisElement E_c(*vec, *c);
        for (int face = 0; face < c.Faces(); ++face) {
          if (E_c.Bnd(face) != -1) continue;
          Point p = GlobalPointOnFace(c, face, RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face);
          Point p_n = neighbour.GlobalToLocal(p);
          ArgyrisElement E_n(*vec, *neighbour);

          // shape functions corresponding to left face corner
          Point np = GlobalPointOnFace(c, face, 0.0);
          int i_c = E_c.NodalPointIndex(np);
          int i_n = E_n.NodalPointIndex(np);
          for (int k = 0; k < 6; ++k) {
            VectorField D_c = E_c.Derivative(p_c, 6 * i_c + k);
            VectorField D_n = E_n.Derivative(p_n, 6 * i_n + k);
            EXPECT_NEAR(D_c[0], D_n[0], testTol);
            EXPECT_NEAR(D_c[1], D_n[1], testTol);
          }

          // shape functions corresponding to right face corner
          np = GlobalPointOnFace(c, face, 1.0);
          i_c = E_c.NodalPointIndex(np);
          i_n = E_n.NodalPointIndex(np);
          for (int k = 0; k < 6; ++k) {
            VectorField D_c = E_c.Derivative(p_c, 6 * i_c + k);
            VectorField D_n = E_n.Derivative(p_n, 6 * i_n + k);
            EXPECT_NEAR(D_c[0], D_n[0], testTol);
            EXPECT_NEAR(D_c[1], D_n[1], testTol);
          }

          // shape functions corresponding to face middle
          np = GlobalPointOnFace(c, face, 0.5);
          i_c = E_c.NodalPointIndex(np);
          i_n = E_n.NodalPointIndex(np);
          VectorField D_c = E_c.Derivative(p_c, 18 + i_c);
          VectorField D_n = E_n.Derivative(p_n, 18 + i_n);
          EXPECT_NEAR(D_c[0], D_n[0], testTol);
          EXPECT_NEAR(D_c[1], D_n[1], testTol);
        }
      }
    }
  }

  void checkNormal() {
    const Mesh &mesh = meshes->fine();
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      ArgyrisElement E_c(*vec, *c);
      for (int face_c = 0; face_c < c.Faces(); ++face_c) {
        if (E_c.Bnd(face_c) != -1) continue;
        cell neighbour = vec->find_neighbour_cell(c, face_c);
        ArgyrisElement E_n(*vec, *neighbour);
        int face_n = (*neighbour).Face(c.Face(face_c));

        EXPECT_EQ(E_c.OrientedNormal(face_c), E_n.OrientedNormal(face_n));
      }
    }
  }
};

TEST_P(TestArgyrisElement, TestContinuity) { checkContinuity(); }

TEST_P(TestArgyrisElement, TestDerivativeContinuity) { checkDerivativeContinuity(); }

TEST_P(TestArgyrisElement, TestNormalConsistency) { checkNormal(); }

TEST_P(TestArgyrisElement, TestMaxK) {
  ArgyrisElement element(*vec, *vec->cells());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(element.get_maxk(i), 6);
  for (int i = 3; i < 6; ++i)
    EXPECT_EQ(element.get_maxk(i), 1);
}

TEST_P(TestArgyrisElement, TestIndexingK) {
  ArgyrisElement element(*vec, *vec->cells());
  for (int m = 0; m < 3; ++m)
    for (int i = 0; i < element.size(); ++i)
      for (int k = 0; k < element.get_maxk(i); ++k)
        EXPECT_EQ(element.indexing_k(i, k, m), element.get_maxk(i) * m + k);
}

TEST_P(TestArgyrisElement, TestIndexingShape) { checkIndexingShape(); }

TEST_P(TestArgyrisElement, TestNodalPointValues) { checkNodalPointValues(); }

INSTANTIATE_TEST_SUITE_P(TestArgyrisElement, TestArgyrisElement,
                         Values("Triangle", "Square2Triangles"));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
