#ifndef TESTDOF_HPP
#define TESTDOF_HPP

#include "IDoF.hpp"
#include "TestEnvironment.hpp"
#include "cells/Cell.hpp"

#include <utility>

struct DoFData {
  int degree;
  vector<Point> nodalPoints;
  vector<short> nodalDoFs;
  vector<vector<int>> nodalPointOnFace;
};

class DoFTest : public TestWithParam<DoFData> {
protected:
  Cell *cell = nullptr;
  IDoF *doF = nullptr;
  std::string dofName = "";

  DoFTest(CELLTYPE type, std::string dofName) : dofName(dofName) {
    cell = CreateCell(type, 0, ReferenceCell(type)->AsVector());
  }

  ~DoFTest() {
    if (doF) delete doF;
    if (cell) delete cell;
  }

  virtual void checkDoFName() { EXPECT_STREQ(doF->Name().c_str(), dofName.c_str()); }

  virtual void checkNumNodalPoints() {
    EXPECT_EQ(GetParam().nodalPoints.size(), doF->NumberOfNodalPoints(*cell));
  }

  virtual void checkNodalPoints() {
    std::vector<Point> nodalPoints = doF->GetNodalPoints(*cell);
    for (int i = 0; i < nodalPoints.size(); ++i)
      EXPECT_POINT_NEAR(GetParam().nodalPoints[i], nodalPoints[i], tolerance);
  }

  virtual void checkNodalDoFs() {
    std::vector<short> nodalDoFs = doF->DoFSizesAtNodalPoints(*cell);
    for (int i = 0; i < nodalDoFs.size(); ++i)
      EXPECT_EQ(GetParam().nodalDoFs[i], nodalDoFs[i]);
  }

  virtual void checkNodalPointsOnFace() {
    for (int face = 0; face < cell->Faces(); ++face)
      EXPECT_EQ(GetParam().nodalPointOnFace[face].size(),
                doF->NumberOfNodalPointsOnFace(*cell, face));
  }

  virtual void checkNodalPointOnFace() {
    for (int face = 0; face < cell->Faces(); ++face)
      for (int k = 0; k < doF->NumberOfNodalPointsOnFace(*cell, face); ++k)
        EXPECT_EQ(GetParam().nodalPointOnFace[face][k],
                  doF->IdOfStoragePointOnFace(*cell, face, k));
  }
};

/// To avoid the same code in each test
#define DOF_TESTS(dofTestClass)                                                                    \
                                                                                                   \
  TEST_P(dofTestClass, DoFNameTest) { checkDoFName(); }                                            \
                                                                                                   \
  TEST_P(dofTestClass, NumNodalPointsTest) { checkNumNodalPoints(); }                              \
                                                                                                   \
  TEST_P(dofTestClass, NodalPointsTest) { checkNodalPoints(); }                                    \
                                                                                                   \
  TEST_P(dofTestClass, NodalDoFsTest) { checkNodalDoFs(); }                                        \
                                                                                                   \
  TEST_P(dofTestClass, NodalPointsOnFaceTest) { checkNodalPointsOnFace(); }                        \
                                                                                                   \
  TEST_P(dofTestClass, NodalPointOnFaceTest) { checkNodalPointOnFace(); }

#endif // TESTDOF_HPP