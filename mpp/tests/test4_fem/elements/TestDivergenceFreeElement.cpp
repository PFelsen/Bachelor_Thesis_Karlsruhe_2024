#include "DivergenceFreeDiscretization.hpp"
#include "DivergenceFreeElement.hpp"
#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"
#include "Vector.hpp"

const double testTol = 1e-12;

#if SpaceDimension >= 2

class TestDivergenceFreeElement : public TestWithParam<std::string> {
  Point GlobalPointOnFace(const cell &c, int face, double t) {
    return (*c).FaceLocalToGlobal(face, Point(t, 0.0));
  }
protected:
  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const DivergenceFreeDiscretization> disc;
  Vector *vec = nullptr;

  TestDivergenceFreeElement() {
    meshes = MeshesCreator(GetParam()).WithPLevel(2).WithLevel(2).CreateShared();
    disc = std::make_shared<const DivergenceFreeDiscretization>(*meshes);
    vec = new Vector(0.0, disc);
  }

  ~TestDivergenceFreeElement() { delete vec; }

  // Placed here to exploit friendship
  void checkIndexingShape() {
    DivergenceFreeElement element(*vec, *vec->cells());
    for (int i = 0; i < element.size(); ++i)
      for (int k = 0; k < element.get_maxk(i); ++k) {
        if (i < 3) {
          EXPECT_EQ(element.indexingShape(i, k), i * 6 + k);
        } else {
          EXPECT_EQ(element.indexingShape(i, k), 15 + i);
        }
      }
  }

  void checkContinuity() {
    for (int cnt = 0; cnt < 10; ++cnt) {
      const Mesh &mesh = meshes->fine();
      for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
        DivergenceFreeElement E_c(*vec, *c);
        for (int face = 0; face < c.Faces(); ++face) {
          if (E_c.Bnd(face) != -1) continue;
          Point p = GlobalPointOnFace(c, face, RandomDouble(0.0, 1.0));
          Point p_c = c.GlobalToLocal(p);
          cell neighbour = vec->find_neighbour_cell(c, face);
          Point p_n = neighbour.GlobalToLocal(p);
          DivergenceFreeElement E_n(*vec, *neighbour);

          // shape functions corresponding to left face corner
          Point np = GlobalPointOnFace(c, face, 0.0);
          int i_c = E_c.NodalPointIndex(np);
          int i_n = E_n.NodalPointIndex(np);
          for (int k = 0; k < 6; ++k) {
            Velocity D_c = E_c.VelocityField(p_c, 6 * i_c + k);
            Velocity D_n = E_n.VelocityField(p_n, 6 * i_n + k);
            EXPECT_NEAR(D_c[0], D_n[0], testTol);
            EXPECT_NEAR(D_c[1], D_n[1], testTol);
          }

          // shape functions corresponding to right face corner
          np = GlobalPointOnFace(c, face, 1.0);
          i_c = E_c.NodalPointIndex(np);
          i_n = E_n.NodalPointIndex(np);
          for (int k = 0; k < 6; ++k) {
            Velocity D_c = E_c.VelocityField(p_c, 6 * i_c + k);
            Velocity D_n = E_n.VelocityField(p_n, 6 * i_n + k);
            EXPECT_NEAR(D_c[0], D_n[0], testTol);
            EXPECT_NEAR(D_c[1], D_n[1], testTol);
          }

          // shape functions corresponding to face middle
          np = GlobalPointOnFace(c, face, 0.5);
          i_c = E_c.NodalPointIndex(np);
          i_n = E_n.NodalPointIndex(np);
          Velocity D_c = E_c.VelocityField(p_c, 18 + i_c);
          Velocity D_n = E_n.VelocityField(p_n, 18 + i_n);
          EXPECT_NEAR(D_c[0], D_n[0], testTol);
          EXPECT_NEAR(D_c[1], D_n[1], testTol);
        }
      }
    }
  }

  void checkNormal() {
    const Mesh &mesh = meshes->fine();
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      DivergenceFreeElement E_c(*vec, *c);
      for (int face_c = 0; face_c < c.Faces(); ++face_c) {
        if (E_c.Bnd(face_c) != -1) continue;
        cell neighbour = vec->find_neighbour_cell(c, face_c);
        ArgyrisElement E_n(*vec, *neighbour);
        int face_n = (*neighbour).Face(c.Face(face_c));

        EXPECT_EQ(E_c.OrientedNormal(face_c), E_n.OrientedNormal(face_n));
      }
    }
  }

  void checkDivergence() {
    const Mesh &mesh = meshes->fine();
    for (cell c = vec->cells(); c != vec->cells_end(); ++c) {
      DivergenceFreeElement element(*vec, *c);
      for (int i = 0; i < element.size(); ++i) {
        for (int k = 0; k < element.get_maxk(i); ++k) {
          for (int q = 0; q < element.nQ(); ++q) {
            Tensor D = element.VelocityFieldGradient(q, i, k);
            EXPECT_NEAR(D[0][0] + D[1][1], 0.0, testTol);
          }
        }
      }
    }
  }
};

TEST_P(TestDivergenceFreeElement, TestContinuity) { checkContinuity(); }

TEST_P(TestDivergenceFreeElement, TestNormalConsistency) { checkNormal(); }

TEST_P(TestDivergenceFreeElement, TestMaxK) {
  DivergenceFreeElement element(*vec, *vec->cells());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(element.get_maxk(i), 6);
  for (int i = 3; i < 6; ++i)
    EXPECT_EQ(element.get_maxk(i), 1);
}

TEST_P(TestDivergenceFreeElement, TestIndexingK) {
  ArgyrisElement element(*vec, *vec->cells());
  for (int m = 0; m < 3; ++m)
    for (int i = 0; i < element.size(); ++i)
      for (int k = 0; k < element.get_maxk(i); ++k)
        EXPECT_EQ(element.indexing_k(i, k, m), element.get_maxk(i) * m + k);
}

TEST_P(TestDivergenceFreeElement, TestIndexingShape) { checkIndexingShape(); }

TEST_P(TestDivergenceFreeElement, TestDivergenceZero) { checkDivergence(); }

INSTANTIATE_TEST_CASE_P(TestDivergenceFreeElement, TestDivergenceFreeElement,
                        Values("Triangle", "Square2Triangles"));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
