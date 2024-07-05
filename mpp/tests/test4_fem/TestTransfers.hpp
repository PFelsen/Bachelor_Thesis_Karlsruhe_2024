#ifndef TESTTRANSFERS_HPP
#define TESTTRANSFERS_HPP

#include "TestEnvironment.hpp"
#include "Vector.hpp"
#include "transfer/ITransfer.hpp"

void FillVectorConstant(Vector &vec, double d) { vec = d; }

void ASSERT_CONSTANTVECTOR(const Vector &vec, double d) {
  for (row r = vec.rows(); r != vec.rows_end(); ++r) {
    for (int i = 0; i < r.n(); ++i) {
      EXPECT_NEAR(vec(r, i), d, tolerance);
    }
  }
}

void FillVectorLinear(Vector &vec, double lower, double upper) {
  for (row r = vec.rows(); r != vec.rows_end(); ++r) {
    auto p = r();
    for (int i = 0; i < r.n(); ++i) {
      vec(r, i) = lower + (upper - lower) * p[i];
    }
  }
}

void ASSERT_LINEARVECTOR(const Vector &vec, double lower, double upper) {
  for (row r = vec.rows(); r != vec.rows_end(); ++r) {
    auto p = r();
    for (int i = 0; i < r.n(); ++i) {
      EXPECT_NEAR(vec(r, i), lower + (upper - lower) * p[i], tolerance);
    }
  }
}

void FillVectorQuadratic(Vector &vec, double lower, double upper) {
  for (row r = vec.rows(); r != vec.rows_end(); ++r) {
    auto p = r();
    for (int i = 0; i < r.n(); ++i) {
      auto x = (p[i] - (upper - lower));
      vec(r, i) = lower + x * x;
    }
  }
}

void ASSERT_QUADRATICVECTOR(const Vector &vec, double lower, double upper) {
  for (row r = vec.rows(); r != vec.rows_end(); ++r) {
    auto p = r();
    for (int i = 0; i < r.n(); ++i) {
      auto x = (p[i] - (upper - lower));
      EXPECT_NEAR(vec(r, i), lower + x * x, tolerance);
    }
  }
}

void ASSERT_RESTRICTION(const Vector &coarse, const Vector &fine, RMatrix transfer) {
  transfer.transpose();
  RVector f(fine.nR());
  for (int i = 0; i < fine.nR(); ++i) {
    f[i] = fine(i, 0);
  }
  RVector restricted = transfer * f;

  for (row r = coarse.rows(); r != coarse.rows_end(); ++r) {
    procset p = coarse.find_procset(r());
    double rowNorm = transfer.row(r.Id()).normOnes();
    if (rowNorm > 0) { EXPECT_NEAR(coarse(r, 0), restricted[r.Id()] / rowNorm, 1e-10); }
  }
}

struct TransferParam {
  int coarseLevel;
  int coarseDegree;
  int fineLevel;
  int fineDegree;
  std::string meshName;
};

static std::vector<TransferParam> TransferParameters(std::string geometryName) {
  return std::vector<TransferParam>{{0, 1, 0, 1, geometryName}, {0, 1, 1, 1, geometryName},
                                    {0, 1, 2, 1, geometryName}, {0, 1, 3, 1, geometryName},
                                    {1, 1, 1, 1, geometryName}, {1, 1, 2, 1, geometryName},
                                    {1, 1, 3, 1, geometryName}, {1, 1, 4, 1, geometryName},
                                    {0, 1, 0, 2, geometryName}, {0, 1, 1, 2, geometryName},
                                    {0, 1, 2, 2, geometryName}, {0, 2, 1, 1, geometryName},
                                    {0, 2, 2, 1, geometryName}, {0, 2, 3, 1, geometryName},
                                    {0, 2, 1, 2, geometryName}, {0, 2, 2, 2, geometryName},
                                    {0, 2, 1, 4, geometryName}, {0, 4, 1, 2, geometryName},
                                    {0, 4, 2, 1, geometryName}};
}

#endif // TESTTRANSFERS_HPP
