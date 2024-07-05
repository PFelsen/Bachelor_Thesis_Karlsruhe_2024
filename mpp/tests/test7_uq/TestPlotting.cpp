#include "Sample.hpp"

#include "DGDiscretization.hpp"
#include "DoFDiscretization.hpp"
#include "LagrangeDiscretization.hpp"
#include "RTLagrangeDiscretization.hpp"

#include "MeshesCreator.hpp"
#include "Plotting.hpp"

#include "TestEnvironment.hpp"

#include <filesystem>


namespace fs = std::filesystem;

std::vector<SampleID> testIds =
    {SampleID(0, 0, false),  SampleID(0, 1, true),  SampleID(1, 2, false),  SampleID(1, 3, true),
     SampleID(2, 4, false),  SampleID(2, 5, true),  SampleID(3, 6, false),  SampleID(3, 7, true),
     SampleID(4, 8, false),  SampleID(4, 9, true),  SampleID(5, 10, false), SampleID(5, 11, true),
     SampleID(6, 12, false), SampleID(6, 13, true), SampleID(0, 1, false),  SampleID(0, 2, true),
     SampleID(1, 3, false),  SampleID(1, 4, true),  SampleID(2, 5, false),  SampleID(2, 6, true),
     SampleID(3, 7, false),  SampleID(3, 8, true),  SampleID(4, 9, false),  SampleID(4, 10, true),
     SampleID(5, 11, false), SampleID(5, 12, true), SampleID(6, 13, false), SampleID(6, 14, true)};

class TestPlotting : public Test {
public:
  std::unique_ptr<Meshes> meshes;

  std::vector<std::shared_ptr<IDiscretization>> discs{};

  std::vector<PlotConfig> plotConfigs;

  TestPlotting() {
    fs::remove_all(Config::GetPlotPath());
    std::filesystem::create_directory(Config::GetPlotPath());

    meshes = MeshesCreator("Square")
                 .WithPLevel(testIds.front().cLevel)
                 .WithLevel(testIds.back().fLevel)
                 .CreateUnique();

    discs.push_back(std::make_shared<LagrangeDiscretization>(*meshes, 0));
    discs.push_back(std::make_shared<LagrangeDiscretization>(*meshes, 1));
    discs.push_back(std::make_shared<LagrangeDiscretization>(*meshes, 2));
    discs.push_back(std::make_shared<LagrangeDiscretization>(*meshes, 3));
    discs.push_back(std::make_shared<DGDiscretization>(*meshes, 0));
    discs.push_back(std::make_shared<DGDiscretization>(*meshes, 1));
    discs.push_back(std::make_shared<DGDiscretization>(*meshes, 2));
    discs.push_back(std::make_shared<DGDiscretization>(*meshes, 3));

    plotConfigs.push_back({.parallelPlotting = true});
    plotConfigs.push_back({.parallelPlotting = false});
  }
};

TEST_F(TestPlotting, TestIdsAndDiscs) {
  for (PlotConfig &pConf : plotConfigs) {
    for (auto &disc : discs) {
      for (auto &id : testIds) {
        SampleSolution solution(disc, id, disc->DiscName());
        mpp::plot(solution.IdString(), pConf) << solution.solution.vector << mpp::endp;
        std::string file = Config::GetPlotPath() + solution.IdString();
        if (pConf.parallelPlotting) {
          file += ".pvtu";
        } else {
          file += ".vtu";
        }
        PPM->Barrier(0);
        EXPECT_TRUE(fs::exists(file))
            << DOUT(id) << DOUT(disc->DiscName()) << DOUT(pConf.parallelPlotting);
        // TODO remove plot in order to use other config
        // TODO config will be ignored if plot already exists
        Plotting::Instance().RemovePlot(solution.IdString());
      }
    }
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithScreenLogging().WithPPM()).RUN_ALL_MPP_TESTS();
}
