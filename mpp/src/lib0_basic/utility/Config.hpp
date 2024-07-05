#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <map>
#include <unordered_map>
#include <unordered_set>

#include "Assertion.hpp"
#include "Logging.hpp"
#include "StringUtil.hpp"

class ConfigEntry {
public:
  std::string value;
  std::string origin;

  explicit ConfigEntry(std::string value, std::string origin = "") :
      value(trim_copy(value)), origin(origin) {}
};

class Config {
  int verbose = 0;

  using Map = std::unordered_map<std::string, ConfigEntry>;
  using StringMap = std::map<std::string, std::string>;
  using Pair = std::pair<const std::string, ConfigEntry>;
  using StringPair = std::pair<const std::string, const std::string>;
  using ChainPair = std::pair<const StringPair, const StringPair>;
  using StringPairMap = std::map<const StringPair, const StringPair>;
protected:
  static std::string configFile;
  static std::string confPath;
  static std::string usedConfFile;
  static std::string usedConfPath;
  static StringPairMap chainConfig;
  static bool saveUsedConf;

  static std::string geoPath;
  static std::string logFile;
  static std::string jsonPath;
  static std::string logPath;
  static std::string dataPath;
  static std::string plotPath;
private:
  static Config *singleton;

  static void addChainConfigPairs(Map &map);

  static Map toMap(const StringMap &conf);

  static Map toMap(int *argc, char **argv);

  static std::string stringFromFile(const std::string &path,
                                    std::unordered_set<std::string> &loadedFiles);

  static std::string toString(const Map &confMap);

  static std::string getConfigFileNameWithPath();

  static std::string getUsedConfigFileNameWithPath();
public:
  static Map toMap(const std::string &confString);

  void PrepareStringForUseConfigEntries(const std::string &defaultName, std::string &name);

  static void parseReferenceConfig(Map &map);

  static std::string GetSourcePath();

  static std::string GetProjectPath();

  static std::string GetBuildPath();

  static void SetConfigFileName(const std::string &configFile);

  static void SetConfPath(const std::string &path);

  static std::string GetConfPath();

  static void SetUsedConfPath(const std::string &path);

  static void SetGeoPath(const std::string &path);

  static std::string GetGeoPath();

  static void SetDataPath(const std::string &path);

  static std::string GetDataPath();

  static void SetPlotPath(const std::string &path);

  static std::string GetPlotPath();

  static void SetLogPath(const std::string &path);

  static void SetJsonPath(const std::string &path);

  static std::string GetLogPath();

  static std::string GetJsonPath();

  static std::string GetLogFileNameWithPath();

  static std::string GetJsonFileNameWithPath();

  static void SaveUsedConf(bool save);

  static void AddChainConfig(const StringPair &triggerPair, const StringPair &resultingPair);

  static void Initialize();

  static void Initialize(const std::string &confFile);

  static void Initialize(int *argc, char **argv);

  static void Initialize(std::initializer_list<StringPair> conf);

  static void Initialize(std::initializer_list<Pair> conf);

  static void Initialize(const StringMap &conf);

  static void Initialize(const Map &conf);

  static void Close();

  static void Clear();

  static Config &Instance();

  inline static bool IsInitialized() { return singleton; }

  static bool Exists(const std::string &key);

  template<typename T>
  static bool Get(const std::string &key, T &out, bool mute = false) {
    return Instance().get<T>(key, out, mute);
  }

  template<typename T>
  static T GetWithoutCheck(const std::string &key, bool mute = false) {
    return Instance().getWithoutCheck<T>(key, mute);
  }

  static void PrintInfo();
protected: // to use constructors in tests
  Map configMap{};
  std::unordered_set<std::string> printedKeys{};

  Config() = delete;

  explicit Config(const Map &conf);

  // public:
  bool exists(const std::string &key);

  template<typename T>
  bool get(const std::string &key, T &out, bool mute = false) {
    if (!configMap.contains(key)) {
      if (!printedKeys.contains(key) || (verbose > 4)) {
        if (!mute) {
          if (verbose >= 2) {
            mout << "reading default: " << key
                 << std::string(std::max(0, 40 - (int)key.length()), '.') << " " << out << endl;
            printedKeys.emplace(key);
          }
        }
      }
      return false;
    }

    out = getWithoutCheck<T>(key, mute);
    return true;
  }

  template<typename T>
  T getWithoutCheck(const std::string &key, bool mute = false) {
    ConfigEntry configEntry = configMap.at(key);

    if (!printedKeys.contains(key) || (verbose > 4)) {
      if (!mute) {
        if (verbose >= 2) {
          mout << "reading config:  " << key
               << std::string(std::max(0, 40 - (int)key.length()), '.') << " " << configEntry.value
               << "\n";
          printedKeys.emplace(key);
        }
      }
    }
    try {
      return parse<T>(configEntry.value);
    } catch (std::exception &e) {
      Exit("Reading " + key + ": Failed to convert '" + configEntry.value + "' into the type '"
           + typeid(T).name() + "'.")
    }
  };

  void printInfo() const;

  void clear();
};

#endif // of #ifndef CONFIG_HPP