#ifndef TESTDISTRIBUTION_HPP
#define TESTDISTRIBUTION_HPP

#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"

using namespace std;

struct TestParams {
  std::string meshName;
  std::string distName;
  int plevel = 2;
  int level = 3;
};

class TestDistribution : public TestWithParam<TestParams> {
protected:
  int pLevel{};
  int level{};
  int verbose = 0;

  int maxCommSplit;

  Meshes *meshes;

  explicit TestDistribution() :
      maxCommSplit(PPM->MaxCommSplit()), level(GetParam().level), pLevel(GetParam().plevel),
      meshes(MeshesCreator(GetParam().meshName)
                 .WithDistribute(GetParam().distName)
                 .WithPLevel(pLevel)
                 .WithLevel(level)
                 .Create()) {
    Config::Get("TestDistributionVerbose", verbose);
  }

  void TearDown() override {
    PPM->Barrier(0);
    delete meshes;
  }
};

INSTANTIATE_TEST_SUITE_P(
    TestDistribution, TestDistribution,
    Values(
        TestParams{"Interval", "Stripes"}, TestParams{"Interval", "Stripes_x"},
        TestParams{"Interval", "Stripes_y"}, TestParams{"Interval", "Stripes_z"},
        TestParams{"Interval", "RCB"}, TestParams{"Interval", "RCBx"},
        TestParams{"Interval", "RCBy"}, TestParams{"Interval", "RCBz"},
        TestParams{"Interval", "RCB2D"}, TestParams{"Interval", "RCB2Dx"},
        TestParams{"Interval", "RCB2Dy"}, TestParams{"Interval", "Stripes_Old"},
        TestParams{"Interval", "RCB_Old"}
#if (SpaceDimension >= 2)
        ,
        TestParams{"Triangle", "Stripes"}, TestParams{"Triangle", "Stripes_x"},
        TestParams{"Triangle", "Stripes_y"}, TestParams{"Triangle", "Stripes_z"},
        TestParams{"Triangle", "RCB"}, TestParams{"Triangle", "RCBx"},
        TestParams{"Triangle", "RCBy"}, TestParams{"Triangle", "RCBz"},
        TestParams{"Triangle", "RCB2D"}, TestParams{"Triangle", "RCB2Dx"},
        TestParams{"Triangle", "RCB2Dy"}, TestParams{"Triangle", "Stripes_Old"},
        TestParams{"Triangle", "RCB_Old"}, TestParams{"Square", "Stripes"},
        TestParams{"Square", "Stripes_x"}, TestParams{"Square", "Stripes_y"},
        TestParams{"Square", "Stripes_z"}, TestParams{"Square", "RCB"},
        TestParams{"Square", "RCBx"}, TestParams{"Square", "RCBy"}, TestParams{"Square", "RCBz"},
        TestParams{"Square", "RCB2D"}, TestParams{"Square", "RCB2Dx"},
        TestParams{"Square", "RCB2Dy"}, TestParams{"Square", "RCBG"}, TestParams{"Square", "RCBGx"},
        TestParams{"Square", "RCBGy"}, TestParams{"Square", "Stripes_Old"},
        TestParams{"Square", "RCB_Old"}
#endif
#if (SpaceDimension >= 3)
        ,
        TestParams{"Tetrahedron", "Stripes"}, TestParams{"Tetrahedron", "Stripes_x"},
        TestParams{"Tetrahedron", "Stripes_y"}, TestParams{"Tetrahedron", "Stripes_z"},
        TestParams{"Tetrahedron", "RCB"}, TestParams{"Tetrahedron", "RCBx"},
        TestParams{"Tetrahedron", "RCBy"}, TestParams{"Tetrahedron", "RCBz"},
        TestParams{"Tetrahedron", "RCB2D"}, TestParams{"Tetrahedron", "RCB2Dx"},
        TestParams{"Tetrahedron", "RCB2Dy"}, TestParams{"Tetrahedron", "Stripes_Old"},
        TestParams{"Tetrahedron", "RCB_Old"}, TestParams{"Hexahedron", "Stripes"},
        TestParams{"Hexahedron", "Stripes_x"}, TestParams{"Hexahedron", "Stripes_y"},
        TestParams{"Hexahedron", "Stripes_z"}, TestParams{"Hexahedron", "RCB"},
        TestParams{"Hexahedron", "RCBx"}, TestParams{"Hexahedron", "RCBy"},
        TestParams{"Hexahedron", "RCBz"}, TestParams{"Hexahedron", "RCB2D"},
        TestParams{"Hexahedron", "RCB2Dx"}, TestParams{"Hexahedron", "RCB2Dy"},
        TestParams{"Hexahedron", "RCBG"}, TestParams{"Hexahedron", "RCBGx"},
        TestParams{"Hexahedron", "RCBGy"}, TestParams{"Hexahedron", "RCBGz"},
        TestParams{"Hexahedron", "Stripes_Old"}, TestParams { "Hexahedron", "RCB2" }
#endif
#if (defined USE_SPACETIME && SpaceDimension >= 2)
//        , TODO: does this even make sense for Spacetime. Timelevel is always -1!
//        TestParams{"SpaceTimeSquare", "Stripes",3,4},
//        TestParams{"SpaceTimeSquare", "st-opt",3,4},
//        TestParams{"SpaceTimeSquare", "RCB",3,4},
//        TestParams{"SpaceTimeSquare", "deformed",3,4},
//        TestParams{"SpaceTimeSquare", "deformed_optimized",3,4},
//        TestParams{"SpaceTimeSquare", "optimized_doubleslit",3,4},
//        TestParams{"SpaceTimeSquare", "time_stripes",3,4},
//        TestParams{"SpaceTimeSquare", "x_stripes",3,4},
//        TestParams{"SpaceTimeSquare", "y_stripes",3,4}
#endif
        ));

#endif // TESTDISTRIBUTION_HPP
