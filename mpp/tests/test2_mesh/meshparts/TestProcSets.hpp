#ifndef TESTPROCSET_HPP
#define TESTPROCSET_HPP

#include "LevelPair.hpp"
#include "MeshesCreator.hpp"

#include "TestEnvironment.hpp"

/*
 * Only stripe distributions are tested
 */


class TestProcSets : public TestWithParam<std::string> {
protected:
  std::string distName;

  int commSplit;

  int pLevel = 2;

  int level = 3;

  Meshes *meshes;

  const Mesh &coarseMesh;

  const Mesh &fineMesh;

  TestProcSets(const std::string &meshName, int commSplit = 0) :
      commSplit(commSplit), distName(GetParam()), meshes(MeshesCreator(meshName)
                                                             .WithDistribute(distName)
                                                             .WithPLevel(pLevel)
                                                             .WithLevel(level)
                                                             .Create()),
      coarseMesh((*meshes)[meshes->CoarseLevel().WithCommSplit(commSplit)]),
      fineMesh((*meshes)[meshes->FineLevel().WithCommSplit(commSplit)]) {}

  bool ProcSetSizeLessEqThan(const Mesh &mesh, int lessEqThan) {
    for (procset p = mesh.procsets(); p != mesh.procsets_end(); p++)
      if (p() != Infty && p.size() > lessEqThan) {
        std::cout << "Proc: " << PPM->Proc(commSplit) << " ProcSetSize: " << p.size()
                  << " Should be smaller than " << lessEqThan << endl;
        return false;
      }
    return true;
  }

  void TearDown() override {
    PPM->Barrier(0);
    delete meshes;
  }
};

/*
 * Interval
 */

class TestProcSetsInterval : public TestProcSets {
protected:
  TestProcSetsInterval(int commSplit = 0) : TestProcSets("Interval", commSplit) {}
};

class TestProcSetsIntervalWithSplit : public TestProcSetsInterval {
protected:
  TestProcSetsIntervalWithSplit() : TestProcSetsInterval(1) {}
};

class TestProcSetsIntervalWithDoubleSplit : public TestProcSetsInterval {
protected:
  TestProcSetsIntervalWithDoubleSplit() : TestProcSetsInterval(2) {}
};

#define PROCSETS_INTERVAL_TESTS(ProcSetsIntervalTestClass)                                         \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestProcSets, ProcSetsIntervalTestClass,                                \
                           Values("Stripes_x", "Stripes_y", "Stripes_z"));                         \
                                                                                                   \
  TEST_P(ProcSetsIntervalTestClass, TestProcSetsCount) {                                           \
    EXPECT_EQ(PPM->Size(commSplit) - 1, coarseMesh.ProcSetsCountGeometryWithoutInfty());           \
    EXPECT_EQ(PPM->Size(commSplit) - 1, fineMesh.ProcSetsCountGeometryWithoutInfty());             \
  }                                                                                                \
                                                                                                   \
  TEST_P(ProcSetsIntervalTestClass, TestProcSetsSize) {                                            \
    EXPECT_TRUE(ProcSetSizeLessEqThan(coarseMesh, 2));                                             \
    EXPECT_TRUE(ProcSetSizeLessEqThan(fineMesh, 2));                                               \
  }

/*
 * Square
 */

class TestProcSetsSquare : public TestProcSets {
protected:
  TestProcSetsSquare(int commSplit = 0) : TestProcSets("Square", commSplit) {}
};

class TestProcSetsSquareWithSplit : public TestProcSetsSquare {
protected:
  TestProcSetsSquareWithSplit() : TestProcSetsSquare(1) {}
};

class TestProcSetsSquareWithDoubleSplit : public TestProcSetsSquare {
protected:
  TestProcSetsSquareWithDoubleSplit() : TestProcSetsSquare(2) {}
};

#define PROCSETS_SQUARE_TESTS(ProcSetsSquareTestClass)                                             \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestProcSets, ProcSetsSquareTestClass,                                  \
                           Values("Stripes_x", "Stripes_y", "Stripes_z"));                         \
                                                                                                   \
  TEST_P(ProcSetsSquareTestClass, TestProcSetsCount) {                                             \
    EXPECT_EQ((PPM->Size(commSplit) - 1) * (pow(2, pLevel + 1) + 1),                               \
              coarseMesh.ProcSetsCountGeometryWithoutInfty());                                     \
    EXPECT_EQ((PPM->Size(commSplit) - 1) * (pow(2, level + 1) + 1),                                \
              fineMesh.ProcSetsCountGeometryWithoutInfty());                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(ProcSetsSquareTestClass, TestProcSetsSize) {                                              \
    EXPECT_TRUE(ProcSetSizeLessEqThan(coarseMesh, 2));                                             \
    EXPECT_TRUE(ProcSetSizeLessEqThan(fineMesh, 2));                                               \
  }

/*
 * Hexahedron
 */

class TestProcSetsHexahedron : public TestProcSets {
protected:
  TestProcSetsHexahedron(int commSplit = 0) : TestProcSets("Hexahedron", commSplit) {}

  int ProcSetsCount(int l) {
    int numFaces = pow(pow(2, l), 2);
    int numVertices = pow((pow(2, l) + 1), 2);
    int numEdges = pow(2, l) * (pow(2, l) + 1) * 2;
    return numFaces + numEdges + numVertices;
  }
};

class TestProcSetsHexahedronWithSplit : public TestProcSetsHexahedron {
protected:
  TestProcSetsHexahedronWithSplit() : TestProcSetsHexahedron(1) {}
};

class TestProcSetsHexahedronWithDoubleSplit : public TestProcSetsHexahedron {
protected:
  TestProcSetsHexahedronWithDoubleSplit() : TestProcSetsHexahedron(2) {}
};

#define PROCSETS_HEXAHEDRON_TESTS(DistributeHexahedronTestClass)                                   \
                                                                                                   \
  INSTANTIATE_TEST_SUITE_P(TestProcSets, DistributeHexahedronTestClass,                            \
                           Values("Stripes_x", "Stripes_y", "Stripes_z"));                         \
  TEST_P(DistributeHexahedronTestClass, TestProcSetsCount) {                                       \
    EXPECT_EQ((PPM->Size(commSplit) - 1) * ProcSetsCount(pLevel),                                  \
              coarseMesh.ProcSetsCountGeometryWithoutInfty());                                     \
    EXPECT_EQ((PPM->Size(commSplit) - 1) * ProcSetsCount(level),                                   \
              fineMesh.ProcSetsCountGeometryWithoutInfty());                                       \
  }                                                                                                \
                                                                                                   \
  TEST_P(DistributeHexahedronTestClass, TestProcSetsSize) {                                        \
    EXPECT_TRUE(ProcSetSizeLessEqThan(coarseMesh, 2));                                             \
    EXPECT_TRUE(ProcSetSizeLessEqThan(fineMesh, 2));                                               \
  }


#endif // TESTPROCSET_HPP
