#include "Config.hpp"
#include <json.hpp>
#include "M_IOFiles.hpp"
#include "Parallel.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <json.hpp>

void Config::PrepareStringForUseConfigEntries(const std::string &defaultName, std::string &name) {
  if (name.substr(0, 17) == "UseConfigEntries[" && name.back() == ']') {
    int detail = 0;
    get("UseConfigEntriesDetail", detail, true);
    if (name.size() > 18) {
      std::vector<std::string> keys = split(name.substr(17, name.size() - 18), ",");
      if (keys.empty()) {
        WarningOnMaster("No keys found for UseConfigEntries! Using default name:" + defaultName
                        + " String that was parsed: " + name.substr(17, name.size() - 18)) name =
            defaultName;
        return;
      } else {
        name = "";
      }
      for (auto key : keys) {
        std::string value;
        get(key, value, true);
        if (detail > 0) { name.append(key).append("="); }
        name.append(value).append("_");
      }
      name.erase(name.end() - 1);
    } else {
      WarningOnMaster(
          "No keys found for UseConfigEntries it seems to have no arguments! Using default name:"
          + defaultName) name = defaultName;
      return;
    }
  }
}

using namespace std;

std::string Config::configFile = "m++.conf";
std::string Config::confPath = Config::GetProjectPath() + "/conf/";
std::string Config::usedConfFile = "";
std::string Config::usedConfPath = Config::GetBuildPath() + "/usedconf/";
Config::StringPairMap Config::chainConfig{};
bool Config::saveUsedConf = true;

std::string Config::geoPath = Config::GetProjectPath() + "/conf/geo/";
std::string defaultLogFile = "logfile";
std::string Config::logFile = defaultLogFile;
std::string Config::logPath = Config::GetBuildPath() + "/log/";
std::string Config::jsonPath = Config::GetBuildPath() + "/json/";
std::string Config::dataPath = Config::GetBuildPath() + "/data/";
std::string Config::plotPath = Config::GetBuildPath() + "/data/vtu/";

Config *Config::singleton = nullptr;

void addProgramArgumentsToMap(const int *argc, char **argv,
                              std::unordered_map<std::string, ConfigEntry> &configMap) {
  for (int i = 1; i < *argc; i++) {
    char *line_ptr = *(argv + i);
    string line(line_ptr);
    int equalPos = line.find('=');
    if (equalPos == string::npos) {
      cout << "Skipping entry \'" << line << "\' from command-line: no equal sign found." << endl;
      continue;
    }
    string key = line.substr(0, equalPos);
    trim(key);
    string value = line.substr(equalPos + 1);
    trim(value);

    if (configMap.find(key) != configMap.end()) {
      // Overwriting already inserted key
      configMap.at(key) = ConfigEntry(value);
    } else {
      // Creating new entry with key
      configMap.emplace(make_pair(key, ConfigEntry(value)));
    }
  }
}

void Config::addChainConfigPairs(Map &configMap) {
  for (const std::pair<const StringPair, const StringPair> &chainConfPair : Config::chainConfig) {
    const StringPair &triggerPair = chainConfPair.first;
    /*cout << chainConfPair.first.first << " "
              << chainConfPair.first.second << " -> "
              << chainConfPair.second.first << " "
              << chainConfPair.second.second << endl;*/

    if (configMap.find(triggerPair.first) != configMap.end()
        && configMap.at(triggerPair.first).value == triggerPair.second) {
      const StringPair &resultingPair = chainConfPair.second;

      string origin = "ChainConfig from Pair ";
      origin += "(" + triggerPair.first + "," + triggerPair.second + ")";

      ConfigEntry value(resultingPair.second, origin);

      string key = resultingPair.first;

      if (configMap.find(key) != configMap.end()) {
        // Overwriting already inserted key
        configMap.at(key) = ConfigEntry(value);
      } else {
        // Creating new entry with key
        configMap.emplace(make_pair(key, ConfigEntry(value)));
      }
    }
  }
}

Config::Map Config::toMap(const StringMap &conf) {
  Map confMap;
  for (const StringPair &pair : conf) {
    confMap.emplace(pair.first, ConfigEntry(pair.second));
  }
  return confMap;
}

Config::Map Config::toMap(const std::string &confString) {
  // content is in format ((g1)=(g2);)*
  Map confMap;
  for (const string &line : split(confString, ";")) {
    if (line.empty()) { continue; }
    size_t equalPos = line.find('=');
    if (equalPos == string::npos) {
      cout << ">" << line << "<" << endl;
      cout << "No equal sign in a line of the settings file! Skipping entry..." << endl;
      continue;
    }
    string key = line.substr(0, equalPos);
    trim(key);
    string value = line.substr(equalPos + 1);
    trim(value);
    confMap.emplace(key, ConfigEntry(value));
  }
  return confMap;
}

Config::Map Config::toMap(int *argc, char **argv) {
  if (ParallelProgrammingModel::IsInitialized() && !PPM->Master(0)) {
    // When running on multiple processes, non-master processes are informed
    // about the config to be used by master.
    int len;
    PPM->BCastOnCommSplit(len, 0);
    // Note: content is NOT a zero-terminated string.
    char content[len];
    PPM->Broadcast(content, len, 0);
    return toMap(string(content, len));
  }


  // Remove all leading arguments that look like config files,
  // treating the last one as the config file that should actually be used
  // FIXME: Simplify this!
  bool isJsonConfig = false;
  while (*argc > 1 && !contains(argv[1], "=")) {
    // Configuration parameters can either be read from .conf files or
    // from a JSON log produced during a previous run.
    // These are distinguished by their file extension:
    // * JSON logs end with ".json"
    // * Configuration files are *always* specified without an extension
    std::string configFileName = std::string(argv[1]);
    isJsonConfig = configFileName.ends_with(".json");
    if (!isJsonConfig) { configFileName += ".conf"; }

    SetConfigFileName(configFileName);

    for (int i = 1; i < *argc - 1; i++) {
      argv[i] = argv[i + 1];
    }
    --(*argc);
  }

  Map configMap;
  if (isJsonConfig) {
    // Load the json log file
    std::string path = getConfigFileNameWithPath();
    std::ifstream ifs(path);
    if (!ifs.good()) Exit("Failed to open " + path);
    nlohmann::json logFile = nlohmann::json::parse(ifs);

    if (!logFile.contains("Config Info"))
      Exit(path
           + " does not contain a \"Config Info\" entry, did you forget to call "
             "Config::PrintInfo() ?");

    // Convert the json object to a regular configuration map
    configMap.reserve(logFile["Config Info"].size());
    for (auto &el : logFile["Config Info"].items()) {
      if (el.key() != "time")
        configMap.emplace(std::make_pair(el.key(), el.value().get<std::string>()));
    }
  } else {
    // Keep track of the configuration files that were loaded.
    // Config files can reference each other so in theory, circular references
    // are possible and need to be detected.
    std::unordered_set<std::string> loadedFiles;

    // Load the config file, resolving "loadconf" references if present
    std::string content = stringFromFile(getConfigFileNameWithPath(), loadedFiles);

    // Parse the configuration file
    configMap = toMap(content);
  }

  // Combine the config with parameters set using commandline arguments
  addProgramArgumentsToMap(argc, argv, configMap);

  addChainConfigPairs(configMap);

  parseReferenceConfig(configMap);

  if (ParallelProgrammingModel::IsInitialized()) {
    // If we are running on multiple processes, then we need to tell everyone else
    // which configuration values they should be using.

    // FIXME: We can send strings direcly through the buffer, we likely don't need
    // to treat them as raw bytes here
    string cont = toString(configMap);
    int len = cont.size();
    const char *content = cont.c_str();
    char c[len + 1];
    strcpy(c, content);
    PPM->BCastOnCommSplit(len, 0);
    PPM->Broadcast(c, len, 0);
  }

  return configMap;
}

// FIXME: This is parsing the config file, checking for "loadconf" entries, resolving them if
// necessary and then serializing the entire thing back to a single string so we can parse it
// *again*. This is inefficient.
string Config::stringFromFile(const string &path, unordered_set<string> &loadedFiles) {
  // loadedFiles contains all previously loaded conf-files
  // to prevent endless loading if cycles are present

  std::ifstream file(path);
  if (!file.good()) { Exit("No text-file found at: " + path) }
  std::string content;
  std::string line;

  while (getline(file, line)) {
    // check for comments
    string commentDelimiters[] = {"#", "//", ";"};
    for (const string &commentDelimeter : commentDelimiters) {
      int commentPos = line.find(commentDelimeter);
      if (commentPos != string::npos) {
        line = line.substr(0, commentPos);
        trim(line);
      }
    }
    if (line.length() == 0) { continue; }

    size_t equalPos = line.find('=');
    if (equalPos == string::npos) { continue; }
    string key = line.substr(0, equalPos);
    trim(key);
    string value = line.substr(equalPos + 1);
    trim(value);

    if (key == "loadconf") {
      // If file was read previously, skip parsing it again.
      if (loadedFiles.find(value) != loadedFiles.end()) { continue; }
      content += Config::stringFromFile(Config::GetConfPath() + value, loadedFiles);
    } else {
      content += key.append("=").append(value).append(";");
    }
  }
  file.close();
  return content;
}

std::string Config::toString(const Map &confMap) {
  std::string content;
  for (const auto &p : confMap) {
    content += p.first + "=" + p.second.value + ";";
  }
  return content;
}

std::string Config::getConfigFileNameWithPath() { return confPath + configFile; }

std::string Config::getUsedConfigFileNameWithPath() { return usedConfPath + usedConfFile; }

void Config::parseReferenceConfig(Map &configMap) {
  for (auto &[key, ce] : configMap) {
    auto &value = ce.value;
    if (value[0] == '&') {
      auto iter = configMap.find(value.substr(1));
      if (iter != configMap.end()) {
        ce = iter->second;
      } else {
        std::cerr << "Reference-key " << value.substr(1) << " for key " << key << " not found!"
                  << std::endl;
      }
    }
  }
}

std::string Config::GetSourcePath() { return std::string(ProjectSourceDir); }

std::string Config::GetProjectPath() {
  if (std::string(ProjectDir).empty()) { return GetSourcePath(); }
  return GetSourcePath() + "/" + std::string(ProjectDir);
}

std::string Config::GetBuildPath() { return std::string(ProjectBuildDir); }

void Config::SetConfigFileName(const std::string &configFileName) {
  if (singleton) THROW("Config already initialized! Call Config::SetConfigFileName in advance")
  configFile = configFileName;
}

void Config::SetDataPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetDataPath in advance")
  dataPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(dataPath)) { std::filesystem::create_directory(dataPath); }
}

string Config::GetDataPath() { return dataPath; }

void Config::SetConfPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetConfPath in advance")
  confPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(confPath)) { std::filesystem::create_directory(confPath); }
}

string Config::GetConfPath() { return confPath; }

void Config::SetUsedConfPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetUsedConfPath in advance")
  usedConfPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(usedConfPath)) {
    std::filesystem::create_directory(usedConfPath);
  }
}

void Config::SetGeoPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetGeoPath in advance")
  geoPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(geoPath)) { std::filesystem::create_directory(geoPath); }
}

string Config::GetGeoPath() { return geoPath; }

void Config::SetPlotPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetPlotPath in advance")
  plotPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(plotPath)) { std::filesystem::create_directory(plotPath); }
}

string Config::GetPlotPath() { return plotPath; }

void Config::SetLogPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetLogPath in advance")
  logPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(logPath)) { std::filesystem::create_directory(logPath); }
}

void Config::SetJsonPath(const std::string &path) {
  if (singleton) THROW("Config already initialized! Call Config::SetJsonPath in advance")
  jsonPath = makeStringPathLike(path);
  if (!std::filesystem::is_directory(jsonPath)) { std::filesystem::create_directory(jsonPath); }
}

string Config::GetLogPath() { return logPath; }

string Config::GetJsonPath() { return jsonPath; }

std::string Config::GetLogFileNameWithPath() { return logPath + logFile; }

std::string Config::GetJsonFileNameWithPath() { return jsonPath + logFile; }

void Config::SaveUsedConf(bool save) {
  if (singleton) THROW("Config already initialized! Call Config::SaveUsedConf in advance")
  saveUsedConf = save;
}

void Config::AddChainConfig(const StringPair &triggerPair, const StringPair &resultingPair) {
  chainConfig.emplace(triggerPair, resultingPair);
}

void Config::Initialize() { Initialize(Map{}); }

void Config::Initialize(const std::string &confFile) {
  Config::SetConfigFileName(confFile);
  int dummy = 0;
  Initialize(&dummy, nullptr);
}

void Config::Initialize(int *argc, char **argv) { Initialize(toMap(argc, argv)); }

void Config::Initialize(std::initializer_list<StringPair> conf) {
  Initialize(StringMap(conf.begin(), conf.end()));
}

void Config::Initialize(std::initializer_list<Pair> conf) {
  Initialize(Map(conf.begin(), conf.end()));
}

void Config::Initialize(const Config::StringMap &conf) { Initialize(toMap(conf)); }

void Config::Initialize(const Config::Map &conf) {
  if (singleton) THROW("Config already initialized!")
  singleton = new Config(conf);
}

void Config::Close() {
  if (singleton) delete singleton;
  singleton = nullptr;
  chainConfig.clear();
  logFile = defaultLogFile;
}

void Config::Clear() { Instance().clear(); }

Config &Config::Instance() {
  if (!IsInitialized()) THROW("Config not initialized!")
  return *singleton;
}

bool Config::Exists(const std::string &key) { return Instance().exists(key); }

void Config::PrintInfo() { Instance().printInfo(); }

Config::Config(const Map &conf) : configMap(conf) {
  get("ConfigVerbose", verbose, true);

  string path = "";
  get("GeoPath", path, true);
  if (path != "") SetGeoPath(path);
  path = "";
  get("LogPath", path, true);
  if (path != "") SetLogPath(path);
  path = "";
  get("JsonPath", path, true);
  if (path != "") SetJsonPath(path);
  path = "";
  get("DataPath", path, true);
  if (path != "") SetDataPath(path);
  path = "";
  get("PlotPath", path, true);
  if (path != "") SetPlotPath(path);
  path = "";
  get("UsedConfPath", path, true);
  if (path != "") SetUsedConfPath(path);

  get("logfile", logFile, true);
  PrepareStringForUseConfigEntries(defaultLogFile, logFile);
  usedConfFile = logFile;
  usedConfFile.append(".used.conf");

  if (saveUsedConf && (!ParallelProgrammingModel::IsInitialized() || PPM->Master(0))) {
    M_ofstream out(getUsedConfigFileNameWithPath().c_str(), "used.conf");
    out << "#------------------------------------------------------------------------------" << endl
        << "# Used paths" << endl
        << "#------------------------------------------------------------------------------" << endl
        << "# SourcePath=  " << GetSourcePath() << endl
        << "# ProjectPath= " << GetProjectPath() << endl
        << "# ConfPath=    " << GetConfPath() << endl
        << "# GeoPath=     " << GetGeoPath() << endl
        << "# BuildPath=   " << GetBuildPath() << endl
        << "# LogPath=     " << GetLogPath() << endl
        << "# JsonPath=     " << GetJsonPath() << endl
        << "# DataPath=    " << GetDataPath() << endl
        << "# PlotPath=    " << GetPlotPath() << endl
        << "# UsedConfPath=" << usedConfPath << endl
        << "#------------------------------------------------------------------------------" << endl
        << endl;
    for (const Pair p : configMap) {
      out << p.first << "=" << p.second.value << endl;
    }
    out.close();
  }
}

bool Config::exists(const std::string &key) { return configMap.find(key) != configMap.end(); }

// Logs the specified configuration values in a way that they can later be used to start new
// computations.
void Config::printInfo() const {
  if (verbose < 1) return;

  vector<PrintInfoEntry<std::string>> entries;
  entries.reserve(configMap.size());

  for (const auto &p : configMap)
    entries.push_back(PrintInfoEntry<std::string>(p.first, p.second.value, 1));
  std::sort(entries.begin(), entries.end());

  mout.PrintInfo("Config", verbose, entries);
}

void Config::clear() { printedKeys.clear(); }
