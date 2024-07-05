#include "CoarseGeometry.hpp"
#include "Config.hpp"
#include "m++.hpp"

#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  if (argc < 2) {
    Warning("No path was provided. Call ConvertGeoFiles with proper file path") return 0;
  }
  std::string geoFolder = argv[1];
  std::string confName = "default";
  strcpy(argv[1], confName.c_str());

  Config::SetConfPath(std::string(ProjectMppDir) + "/tools/geometries/");
  Config::SetGeoPath(geoFolder);
  Mpp::initialize(&argc, argv);
  mout.setFileEnabled(false);

  std::cout << geoFolder << std::endl;
  for (const auto &file : std::filesystem::directory_iterator(geoFolder)) {
    // get filename
    std::string base_filename = file.path().filename();

    // remove extension from filename
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string filename = base_filename.substr(0, p);
    std::cout << filename << std::endl;
    CoarseGeometry geo(filename);
  }

  return 0;
}