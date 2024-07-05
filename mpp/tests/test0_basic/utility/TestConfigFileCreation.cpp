#include <filesystem>
#include "TestEnvironment.hpp"
#include "utility/Config.hpp"

class ConfigPathCreation : public TestWithParam<std::string> {
  std::string entryname = "";
  std::string pathname = "";
public:
  ConfigPathCreation() : entryname(GetParam()), pathname(GetParam() + "Path") {}

  void SetUp() override {}

  void TearDown() override { Config::Close(); }

  void CheckIfPathExists() {
    std::map<std::string, std::string> configMap{{entryname, pathname}};
    Config::Initialize(configMap);
    EXPECT_EQ(std::filesystem::is_directory(pathname), true);
    std::filesystem::remove(pathname);
  }
};

TEST_P(ConfigPathCreation, TestPathCreation) { CheckIfPathExists(); }

INSTANTIATE_TEST_SUITE_P(ConfigPathCreationTest, ConfigPathCreation,
                         Values("LogPath", "GeoPath", "DataPath", "PlotPath"));
