#include "TestTransportEquality_ST_TS.hpp"

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfPath(std::string(ProjectMppDir) + "/conf/")
                     .WithGeoPath(string(ProjectMppDir) + "/conf/geo/")
                     .
                 // WithoutDefaultConfig().
                 WithConfigEntry("T", 1)
                     .WithConfigEntry("dt", 0.01)
                     .WithConfigEntry("level", 0)
                     .WithConfigEntry("plevel", 0)
                     .WithConfigEntry("ClearDistribution", 0)
                     .WithConfigEntry("LinearVerbose", 1)
                     .WithScreenLogging()
                     .WithFileLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}