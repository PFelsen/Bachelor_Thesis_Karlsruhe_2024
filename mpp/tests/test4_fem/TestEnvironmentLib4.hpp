#ifndef TESTENVIRONMENTLIB4_HPP
#define TESTENVIRONMENTLIB4_HPP

#include "Vector.hpp"
#include "test3_disc/TestEnvironmentLib3.hpp"

void EXPECT_MPPVECTOR_EQ(const Vector &vec1, const Vector &vec2) {
  EXPECT_EQ(vec1.GetMesh().Name(), vec2.GetMesh().Name());
  EXPECT_EQ(vec1.SpaceLevel(), vec2.SpaceLevel());
  EXPECT_EQ(vec1.TimeLevel(), vec2.TimeLevel());
  EXPECT_EQ(typeid(vec1.GetDisc()), typeid(vec2.GetDisc()));

  for (row r = vec1.rows(); r != vec1.rows_end(); ++r) {
    for (int i = 0; i < r.size(); ++i) {
      EXPECT_DOUBLE_EQ(vec1(r, i), vec2(r, i));
    }
  }
}

void ASSERT_MPPVECTOR_EQ(const Vector &vec1, const Vector &vec2) {
  ASSERT_EQ(vec1.GetMesh().Name(), vec2.GetMesh().Name());
  ASSERT_EQ(vec1.SpaceLevel(), vec2.SpaceLevel());
  ASSERT_EQ(vec1.TimeLevel(), vec2.TimeLevel());
  ASSERT_EQ(typeid(vec1.GetDisc()), typeid(vec2.GetDisc()));

  for (row r = vec1.rows(); r != vec1.rows_end(); ++r) {
    for (int i = 0; i < r.n(); ++i) {
      ASSERT_DOUBLE_EQ(vec1(r, i), vec2(r, i));
    }
  }
}

void EXPECT_MPPVECTOR_NEAR(const Vector &vec1, const Vector &vec2, double tol = 1e-8) {
  EXPECT_EQ(vec1.GetMesh().Name(), vec2.GetMesh().Name());
  EXPECT_EQ(vec1.SpaceLevel(), vec2.SpaceLevel());
  EXPECT_EQ(vec1.TimeLevel(), vec2.TimeLevel());
  EXPECT_EQ(typeid(vec1.GetDisc()), typeid(vec2.GetDisc()));

  for (row r = vec1.rows(); r != vec1.rows_end(); ++r) {
    for (int i = 0; i < r.size(); ++i) {
      EXPECT_NEAR(vec1(r, i), vec2(r, i), tol);
    }
  }
}

void ASSERT_MPPVECTOR_NEAR(const Vector &vec1, const Vector &vec2, double tol = 1e-8) {
  ASSERT_EQ(vec1.GetMesh().Name(), vec2.GetMesh().Name());
  ASSERT_EQ(vec1.SpaceLevel(), vec2.SpaceLevel());
  ASSERT_EQ(vec1.TimeLevel(), vec2.TimeLevel());
  ASSERT_EQ(typeid(vec1.GetDisc()), typeid(vec2.GetDisc()));

  for (row r = vec1.rows(); r != vec1.rows_end(); ++r) {
    for (int i = 0; i < r.n(); ++i) {
      ASSERT_NEAR(vec1(r, i), vec2(r, i), tol);
    }
  }
}

#endif // TESTENVIRONMENTLIB4_HPP
