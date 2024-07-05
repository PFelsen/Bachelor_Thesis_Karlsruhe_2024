#include "ArgyrisDiscretization.hpp"
#include "ArgyrisElement.hpp"
#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"
#include "Vector.hpp"

#include <vector>

using std::vector;

const double testTol = 1e-12;

#if SpaceDimension >= 2

struct TestData {
  std::string mesh;
  int numberOfVectors;
};

class TestArgyrisBernsteinRanges : public TestWithParam<TestData> {
protected:
  Meshes *meshes = nullptr;
  std::shared_ptr<const ArgyrisDiscretization> disc;
  std::shared_ptr<const IAArgyrisDiscretization> IAdisc;
  std::vector<Vector *> vec;
  std::vector<double> lambda;

  TestArgyrisBernsteinRanges() :
      vec(GetParam().numberOfVectors), lambda(GetParam().numberOfVectors) {
    meshes = MeshesCreator(GetParam().mesh).WithPLevel(3).WithLevel(3).Create();
    disc = std::make_shared<const ArgyrisDiscretization>(*meshes, 18, 1);
    IAdisc = std::make_shared<const IAArgyrisDiscretization>(*meshes, 18, 1);
    for (int n = 0; n < vec.size(); ++n) {
      vec[n] = new Vector(0.0, disc);
      for (int i = 0; i < vec[n]->size(); ++i) {
        (*vec[n])[i] = RandomDouble();
      }
      vec[n]->SetIADisc(IAdisc);
      lambda[n] = RandomDouble();
    }
  }

  ~TestArgyrisBernsteinRanges() {
    for (int n = 0; n < vec.size(); ++n)
      delete vec[n];
    delete meshes;
  }

  bool testRange() {
    bool check = true;
    if (vec.size() == 1) {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        IABernsteinRanges BR_c(E_c, *vec[0]);
        IAInterval range;
        BR_c.Range(range);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAInterval u = E_c.Value(q, *vec[0]);
          check = check && (u <= range);
        }
      }
    } else {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        std::vector<IABernsteinData> data{};
        for (int n = 0; n < vec.size(); ++n) {
          data.emplace_back(IABernsteinData{lambda[n], *vec[n]});
        }
        IABernsteinRanges BR_c(E_c, data);
        IAInterval range;
        BR_c.Range(range);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAInterval u;
          for (int n = 0; n < vec.size(); ++n) {
            u += IAInterval(lambda[n]) * E_c.Value(q, *vec[n]);
          }
          check = check && (u <= range);
        }
      }
    }
    return PPM->And(check);
  }

  template<int Z>
  bool testRangeGradient() {
    bool check = true;
    if (vec.size() == 1) {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        IABernsteinRanges BR_c(E_c, *vec[0]);
        IAInterval range[2];
        BR_c.RangeGradient(range[0], range[1]);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAVectorField G = E_c.Derivative(q, *vec[0]);
          check = check && (G[Z] <= range[Z]);
        }
      }
    } else {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        std::vector<IABernsteinData> data{};
        for (int n = 0; n < vec.size(); ++n) {
          data.emplace_back(IABernsteinData{lambda[n], *vec[n]});
        }
        IABernsteinRanges BR_c(E_c, data);
        IAInterval range[2];
        BR_c.RangeGradient(range[0], range[1]);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAVectorField G;
          for (int n = 0; n < vec.size(); ++n) {
            G += IAInterval(lambda[n]) * E_c.Derivative(q, *vec[n]);
          }
          check = check && (G[Z] <= range[Z]);
        }
      }
    }
    return PPM->And(check);
  }

  template<int Z1, int Z2>
  bool testRangeHessian() {
    bool check = true;
    if (vec.size() == 1) {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        IABernsteinRanges BR_c(E_c, *vec[0]);
        IAInterval range[2][2];
        BR_c.RangeHessian(range[0][0], range[0][1], range[1][1]);
        range[1][0] = range[0][1];
        for (int q = 0; q < E_c.nQ(); ++q) {
          IASymTensor H = E_c.Hessian(q, *vec[0]);
          check = check && (H(Z1, Z2) <= range[Z1][Z2]);
        }
      }
    } else {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        std::vector<IABernsteinData> data{};
        for (int n = 0; n < vec.size(); ++n) {
          data.emplace_back(IABernsteinData{lambda[n], *vec[n]});
        }
        IABernsteinRanges BR_c(E_c, data);
        IAInterval range[2][2];
        BR_c.RangeHessian(range[0][0], range[0][1], range[1][1]);
        range[1][0] = range[0][1];
        for (int q = 0; q < E_c.nQ(); ++q) {
          IASymTensor H;
          for (int n = 0; n < vec.size(); ++n) {
            H += IAInterval(lambda[n]) * E_c.Hessian(q, *vec[n]);
          }
          check = check && (H(Z1, Z2) <= range[Z1][Z2]);
        }
      }
    }
    return PPM->And(check);
  }

  bool testAbsMax() {
    bool check = true;
    if (vec.size() == 1) {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        IABernsteinRanges BR_c(E_c, *vec[0]);
        double absMax;
        BR_c.AbsMax(absMax);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAInterval u = E_c.Value(q, *vec[0]);
          check = check && (u <= IAInterval(-absMax, absMax));
        }
      }
    } else {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        std::vector<IABernsteinData> data{};
        for (int n = 0; n < vec.size(); ++n) {
          data.emplace_back(IABernsteinData{lambda[n], *vec[n]});
        }
        IABernsteinRanges BR_c(E_c, data);
        double absMax;
        BR_c.AbsMax(absMax);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAInterval u;
          for (int n = 0; n < vec.size(); ++n) {
            u += IAInterval(lambda[n]) * E_c.Value(q, *vec[n]);
          }
          check = check && (u <= IAInterval(-absMax, absMax));
        }
      }
    }
    return PPM->And(check);
  }

  template<int Z>
  bool testAbsMaxGradient() {
    bool check = true;
    if (vec.size() == 1) {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        IABernsteinRanges BR_c(E_c, *vec[0]);
        double absMax[2];
        BR_c.AbsMaxGradient(absMax[0], absMax[1]);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAVectorField G = E_c.Derivative(q, *vec[0]);
          check = check && (G[Z] <= IAInterval(-absMax[Z], absMax[Z]));
        }
      }
    } else {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        std::vector<IABernsteinData> data{};
        for (int n = 0; n < vec.size(); ++n) {
          data.emplace_back(IABernsteinData{lambda[n], *vec[n]});
        }
        IABernsteinRanges BR_c(E_c, data);
        double absMax[2];
        BR_c.AbsMaxGradient(absMax[0], absMax[1]);
        for (int q = 0; q < E_c.nQ(); ++q) {
          IAVectorField G;
          for (int n = 0; n < vec.size(); ++n) {
            G += IAInterval(lambda[n]) * E_c.Derivative(q, *vec[n]);
          }
          check = check && (G[Z] <= IAInterval(-absMax[Z], absMax[Z]));
        }
      }
    }
    return PPM->And(check);
  }

  template<int Z1, int Z2>
  bool testAbsMaxHessian() {
    bool check = true;
    if (vec.size() == 1) {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        IABernsteinRanges BR_c(E_c, *vec[0]);
        double absMax[2][2];
        BR_c.AbsMaxHessian(absMax[0][0], absMax[0][1], absMax[1][1]);
        absMax[1][0] = absMax[0][1];
        for (int q = 0; q < E_c.nQ(); ++q) {
          IASymTensor H = E_c.Hessian(q, *vec[0]);
          check = check && (H(Z1, Z2) <= IAInterval(-absMax[Z1][Z2], absMax[Z1][Z2]));
        }
      }
    } else {
      for (cell c = vec[0]->cells(); c != vec[0]->cells_end(); ++c) {
        const Cell &C = *c;
        IAArgyrisElement E_c(*vec[0], C);
        std::vector<IABernsteinData> data{};
        for (int n = 0; n < vec.size(); ++n) {
          data.emplace_back(IABernsteinData{lambda[n], *vec[n]});
        }
        IABernsteinRanges BR_c(E_c, data);
        double absMax[2][2];
        BR_c.AbsMaxHessian(absMax[0][0], absMax[0][1], absMax[1][1]);
        absMax[1][0] = absMax[0][1];
        for (int q = 0; q < E_c.nQ(); ++q) {
          IASymTensor H;
          for (int n = 0; n < vec.size(); ++n) {
            H += IAInterval(lambda[n]) * E_c.Hessian(q, *vec[n]);
          }
          check = check && (H(Z1, Z2) <= IAInterval(-absMax[Z1][Z2], absMax[Z1][Z2]));
        }
      }
    }
    return PPM->And(check);
  }
};

TEST_P(TestArgyrisBernsteinRanges, RangeTest) { EXPECT_TRUE(testRange()); }

TEST_P(TestArgyrisBernsteinRanges, RangeGradientXTest) { EXPECT_TRUE(testRangeGradient<0>()); }

TEST_P(TestArgyrisBernsteinRanges, RangeGradientYTest) { EXPECT_TRUE(testRangeGradient<1>()); }

TEST_P(TestArgyrisBernsteinRanges, RangeHessianXXTest) { EXPECT_TRUE((testRangeHessian<0, 0>())); }

TEST_P(TestArgyrisBernsteinRanges, RangeHessianXYTest) { EXPECT_TRUE((testRangeHessian<0, 1>())); }

TEST_P(TestArgyrisBernsteinRanges, RangeHessianYXTest) { EXPECT_TRUE((testRangeHessian<1, 0>())); }

TEST_P(TestArgyrisBernsteinRanges, RangeHessianYYTest) { EXPECT_TRUE((testRangeHessian<1, 1>())); }

TEST_P(TestArgyrisBernsteinRanges, AbsMaxTest) { EXPECT_TRUE(testAbsMax()); }

TEST_P(TestArgyrisBernsteinRanges, AbsMaxGradientXTest) { EXPECT_TRUE(testAbsMaxGradient<0>()); }

TEST_P(TestArgyrisBernsteinRanges, AbsMaxGradientYTest) { EXPECT_TRUE(testAbsMaxGradient<1>()); }

TEST_P(TestArgyrisBernsteinRanges, AbsMaxHessianXXTest) {
  EXPECT_TRUE((testAbsMaxHessian<0, 0>()));
}

TEST_P(TestArgyrisBernsteinRanges, AbsMaxHessianXYTest) {
  EXPECT_TRUE((testAbsMaxHessian<0, 1>()));
}

TEST_P(TestArgyrisBernsteinRanges, AbsMaxHessianYXTest) {
  EXPECT_TRUE((testAbsMaxHessian<1, 0>()));
}

TEST_P(TestArgyrisBernsteinRanges, AbsMaxHessianYYTest) {
  EXPECT_TRUE((testAbsMaxHessian<1, 1>()));
}

INSTANTIATE_TEST_SUITE_P(TestBernsteinRanges, TestArgyrisBernsteinRanges,
                         Values(TestData{"Triangle", 1}, TestData{"Square2Triangles", 1},
                                TestData{"Square4Triangles", 1},
                                TestData{"Square4TrianglesTurned", 1}, TestData{"Triangle", 2},
                                TestData{"Triangle", 3}, TestData{"Triangle", 4}));

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
