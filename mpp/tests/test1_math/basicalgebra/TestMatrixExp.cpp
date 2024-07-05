#include "../../TestEnvironment.hpp"
#include "TestBasicAlgebra.hpp"
#include "basicalgebra/RMatrix.hpp"
#include "gtest/gtest.h"
std::vector<double> ThetaM_ = {
    1.495585217958292e-2, 2.539398330063230e-1, 9.504178996162932e-1,
    2.097847961257068e0,  5.371920351148152e0,
};

class ExpRMatrixTest : public testing::TestWithParam<std::pair<bool, std::pair<int, int>>> {
protected:
public:
  ExpRMatrixTest() { srand(124124); }

  void runTriangular(RMatrix &E, RMatrix &EExp, RVector *Vec) {
    double lambda = RandomDouble();
    double alpha = RandomDouble();
    double mu = RandomDouble();
    E[0][0] = lambda;
    E[0][1] = alpha;
    E[1][1] = mu;
    double eps = 1e-16;
    double norm = E.norm();
    auto degs = GetParam().second;
    E *= (ThetaM_[degs.first] - eps) / norm;
    if (degs.first == 3) { E *= pow(2, degs.second); }
    EExp[0][0] = exp(E[0][0]);
    EExp[1][1] = exp(E[1][1]);
    EExp[0][1] = E[0][1] * (exp(E[0][0]) - exp(E[1][1])) / (E[0][0] - E[1][1]);
    if (Vec) {
      E.Exp(*Vec, GetParam().first);
    } else {
      E.Exp(GetParam().first);
    }
  }

  void runFull(RMatrix &E, RMatrix &EExp, RVector *Vec) {
    double det = 5;

    double eps = 1e-16;
    double s = 0.0;
    while (det > 0) {
      E[0][0] = RandomDouble();
      E[0][1] = RandomDouble();
      E[1][1] = RandomDouble();
      E[1][0] = RandomDouble();
      double norm = E.norm();
      auto degs = GetParam().second;
      E *= (ThetaM_[degs.first] - eps) / norm;
      if (degs.first == 3) { E *= pow(2, degs.second); }
      s = E.Trace() / 2.0;
      double ad = (E[0][0] - s) * (E[1][1] - s);
      double bc = E[1][0] * E[0][1];
      det = ad - bc;
    }
    double q = sqrt(-det);
    EExp = (cosh(q) - s * (sinh(q) / q)) * RMatrix(2).Identity() + (sinh(q) / q) * E;
    EExp *= exp(s);
    if (Vec) {
      E.Exp(*Vec, GetParam().first);
    } else {
      E.Exp(GetParam().first);
    }
  }

  void runThreeTimesThree(RMatrix &E, RMatrix &EExp, RVector *Vec) {
    AntisymRMatrix Antisym(3);
    double gamma = 0;
    for (int n = 0; n < 3; ++n) {
      for (int m = n + 1; m < 3; ++m) {
        double d = RandomDouble();
        Antisym(d, n, m);
      }
    }
    double eps = 1e-16;
    RMatrix A(Antisym);
    double norm = A.norm();
    auto degs = GetParam().second;
    A *= (ThetaM_[degs.first] - eps) / norm;
    if (degs.first == 3) { A *= pow(2, degs.second); }
    gamma = A(0, 1) * A(0, 1) + A(0, 2) * A(0, 2) + A(1, 2) * A(1, 2);
    E = A;
    if (Vec) {
      E.Exp(*Vec, GetParam().first);
    } else {
      E.Exp(GetParam().first);
    }
    EExp.Identity();
    gamma = sqrt(gamma);
    RMatrix sum1((sin(gamma) / gamma) * A);
    EExp += sum1;
    RMatrix sum2(((1 - cos(gamma)) / (gamma * gamma)) * A * A);

    EExp += sum2;
    //        Result -= AntisymExp;
  }
};

TEST_P(ExpRMatrixTest,
       ExpMTestTwoTimesTwoTriangular) { // explicit formula from for triangular matrices
  RMatrix E(0.0, 2);
  RMatrix EExp(0.0, 2);
  runTriangular(E, EExp, nullptr);
  E -= EExp;
  EXPECT_NEAR(E.norm() / EExp.norm(), 0, 1e-13);
}

TEST_P(ExpRMatrixTest,
       ExpMTestTwoTimesTwoTriangularVector) { // explicit formula from for triangular matrices
  RMatrix E(0.0, 2);
  RMatrix EExp(0.0, 2);
  RVector Vec(2);
  Vec[0] = 1.0;
  Vec[1] = 0.0;
  runTriangular(E, EExp, &Vec);
  RVector Sol1 = EExp.col(0);
  double diff = (Vec - Sol1).norm() / Sol1.norm();
  Vec[0] = 0.0;
  Vec[1] = 1.0;
  runTriangular(E, EExp, &Vec);
  RVector Sol2 = EExp.col(1);
  diff += (Vec - Sol2).norm() / Sol2.norm();
  EXPECT_NEAR(diff, 0, 1e-13);
}

TEST_P(ExpRMatrixTest, ExpMTestTwoTimesTwoFull) { // explicit formula for 2x2 Matrix (see wiki)
  RMatrix E(0.0, 2);
  RMatrix EExp(0.0, 2);
  runFull(E, EExp, nullptr);
  E -= EExp;
  EXPECT_NEAR(E.norm() / EExp.norm(), 0, 1e-13);
}

TEST_P(ExpRMatrixTest,
       ExpMTestTwoTimesTwoFullVector) { // explicit formula from for triangular matrices
  RMatrix E(0.0, 2);
  RMatrix EExp(0.0, 2);
  RVector Vec(2);
  Vec[0] = 1.0;
  Vec[1] = 0.0;
  runFull(E, EExp, &Vec);
  RVector Sol1 = EExp.col(0);
  double diff = (Vec - Sol1).norm() / Sol1.norm();
  Vec[0] = 0.0;
  Vec[1] = 1.0;
  runFull(E, EExp, &Vec);
  RVector Sol2 = EExp.col(1);
  diff += (Vec - Sol2).norm() / Sol2.norm();
  EXPECT_NEAR(diff, 0, 1e-13);
}

TEST_P(ExpRMatrixTest, ExpMTestThreeTimesThree) { // explicit formula from "some explicit formulas
                                                  // or the matrix exponential"
  RMatrix E(0.0, 3);
  RMatrix EExp(0.0, 3);
  runThreeTimesThree(E, EExp, nullptr);


  E -= EExp;
  EXPECT_NEAR(E.norm() / EExp.norm(), 0, 1e-13);
}

TEST_P(ExpRMatrixTest, ExpMTestThreeTimesThreeVector) { // explicit formula from "some explicit
                                                        // formulas or the matrix exponential"
  RMatrix E(0.0, 3);
  RMatrix EExp(0.0, 3);
  double diff = 0.0;
  RVector Vec(3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Vec[j] = i == j;
    }
    runThreeTimesThree(E, EExp, &Vec);
    RVector Sol = EExp.col(i);
    diff += ((Vec - Sol).norm() / Sol.norm());
  }
  EXPECT_NEAR(diff, 0, 1e-13);
}

INSTANTIATE_TEST_SUITE_P(ExpRMatrixTestTest, ExpRMatrixTest,
                         Values(std::pair<bool, std::pair<int, int>>({false, {0, 1}}),
                                std::pair<bool, std::pair<int, int>>({false, {1, 1}}),
                                std::pair<bool, std::pair<int, int>>({false, {2, 1}}),
                                std::pair<bool, std::pair<int, int>>({false, {3, 1}}),
                                std::pair<bool, std::pair<int, int>>({false, {3, 1}}),
                                std::pair<bool, std::pair<int, int>>({false, {3, 2}}),
                                std::pair<bool, std::pair<int, int>>({false, {3, 3}}),
                                std::pair<bool, std::pair<int, int>>({false, {3, 4}}),
                                std::pair<bool, std::pair<int, int>>({false, {3, 5}}),
                                std::pair<bool, std::pair<int, int>>({true, {0, 1}}),
                                std::pair<bool, std::pair<int, int>>({true, {1, 1}}),
                                std::pair<bool, std::pair<int, int>>({true, {2, 1}}),
                                std::pair<bool, std::pair<int, int>>({true, {3, 1}}),
                                std::pair<bool, std::pair<int, int>>({true, {3, 1}}),
                                std::pair<bool, std::pair<int, int>>({true, {3, 2}}),
                                std::pair<bool, std::pair<int, int>>({true, {3, 3}}),
                                std::pair<bool, std::pair<int, int>>({true, {3, 4}}),
                                std::pair<bool, std::pair<int, int>>({true, {3, 5}})));

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
