#ifndef TESTMESHESPRINTINFO_HPP
#define TESTMESHESPRINTINFO_HPP

#include <LevelPair.hpp>
#include <mesh/MeshesCreator.hpp>

#include "TestEnvironment.hpp"

struct TestConfig {
  std::string meshName;
  std::string distName;
  int commSplit;
  int pLevel;
  int level;
};

class TestMeshesPrintInfo : public TestWithParam<TestConfig> {
protected:
  std::string meshName;

  std::string distName;

  std::string olapName;

  int commSplit;

  int pLevel;

  int level;

  Meshes *meshes;

  TestMeshesPrintInfo() :
      meshName(GetParam().meshName), distName(GetParam().distName), commSplit(GetParam().commSplit),
      pLevel(GetParam().pLevel), level(GetParam().level) {

    meshes = MeshesCreator(meshName)
                 .WithDistribute(distName)
                 .WithPLevel(pLevel)
                 .WithLevel(level)
                 .Create();

    // Ensure meshes with commsplit != 0 are available
    (*meshes)[LevelPair{pLevel, -1, 0, commSplit}];
  }

  void TearDown() override {
    PPM->Barrier(0);
    delete meshes;
  }
};

INSTANTIATE_TEST_SUITE_P(
    TestMeshesPrintInfo, TestMeshesPrintInfo,
    Values(TestConfig{"Interval", "Stripes", 0, 2, 3}, TestConfig{"Interval", "Stripes", 1, 2, 3},
           TestConfig{"Interval", "Stripes", 0, 2, 3}, TestConfig{"Interval", "Stripes", -1, 2, 3}
#if SpaceDimension >= 2
           ,
           TestConfig{"Square", "Stripes", 0, 2, 3}, TestConfig{"Square", "RCB", 1, 2, 3},
           TestConfig{"Square", "RCB", 0, 2, 3}
#if TimeDimension > 0
           ,
           TestConfig{"SpaceTimeSquare", "Stripes", 0, 1, 2},
           TestConfig{"SpaceTimeSquare", "RCB", 0, 1, 2},
           TestConfig{"SpaceTimeSquare", "RCB", 1, 1, 2}
#endif
#endif
#if SpaceDimension >= 3
           ,
           TestConfig{"Hexahedron", "Stripes", 0, 1, 2}, TestConfig{"Hexahedron", "RCB", 0, 1, 2}
#endif
           ));

#endif // TESTMESHESPRINTINFO_HPP
