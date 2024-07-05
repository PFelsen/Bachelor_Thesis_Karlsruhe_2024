#include <memory>

#include "VtuConverter.hpp"

#include "Mesh.hpp"
#include "m++.hpp"

int main(int argc, char **argv) {

  // Reads first program argument as file location
  if (argc < 2) { Warning("No file was provided. Call ConvertVtu with proper file path") return 0; }

  // Load File
  std::string geometryFile = loadFile(argv[1]);
  std::string confName = "default";
  strcpy(argv[1], confName.c_str());

  // Setting up config and mpp to ensure proper mesh generation
  Logging::DisableFileLogging();
  Config::SetConfPath(std::string(ProjectSourceDir) + "/tools/geometries/");
  Config::SetConfigFileName("default.conf");
  Mpp::initialize(&argc, argv);

  // Add variables to converter.
  VtuConverter converter(geometryFile);
  const auto &M = converter.UseCellArrayAsSubdomain(0).UseArrayForCellData(1).Create();

  // Write Geo File to specific location
  converter.WriteGeoFile(M, std::string(ProjectSourceDir) + "/tools/geometries/geo/");

  return 0;
}