
#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <bitset>
#include <complex>
#include <functional>
#include <iomanip>
#include <list>
#include <map>
#include <span>
#include <string>
#include <vector>

#include "Assertion.hpp"
#include "Concepts.hpp"
#include "LoggingBackends.hpp"
#include "M_IOFiles.hpp"
#include "PrintInfo.hpp"
#include "StringUtil.hpp"
#include "TimeDate.hpp"

#ifdef BUILD_IA
#include "IACInterval.hpp"
#include "IAInterval.hpp"
#endif

enum LoggingBackend {
  Text,
  JSON,

  // Keep this as the last variant!
  Count,
};

using EnabledLoggers = std::bitset<static_cast<int>(LoggingBackend::Count)>;

class MasterLogging {
private:
  static std::unique_ptr<MasterLogging> instance;
  EnabledLoggers enabledLoggers;

  MasterLogging(EnabledLoggers enabledLoggers);
public:
  // Whether or not this is the master process
  bool onMaster;

  std::unique_ptr<JsonLogger> jsonLogger;
  std::unique_ptr<LogTextStream<TextLogger>> textLogger;

  // Singleton stuff
  inline static bool isInitialized() { return (bool)instance; }

  inline static MasterLogging &the() {
    if (!MasterLogging::isInitialized()) init();
    return *instance;
  }

  static void init(EnabledLoggers enabledLoggers = {3}) {
    // If this assertion fails then you likely implicitly initialized
    // the logger by accessing mout before.
    Assert(!MasterLogging::isInitialized());

    instance = std::unique_ptr<MasterLogging>(new MasterLogging(enabledLoggers));
  }

  void setFileEnabled(bool to) {
    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->setFileEnabled(to);
  }

  void setScreenEnabled(bool to) {
    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->setScreenEnabled(to);
  }

  EnabledLoggers getEnabledLoggers() { return this->enabledLoggers; }

  // Enable or disable log backends at runtime
  //
  // Note that while backends are (lazily) initialized if necessary,
  // they are never destructed before MasterLogging itself.
  void setEnabledLoggers(EnabledLoggers enabledLoggers) {
    if (enabledLoggers[LoggingBackend::JSON] && this->jsonLogger == nullptr) {
      this->jsonLogger = std::make_unique<JsonLogger>();
    }
    if (enabledLoggers[LoggingBackend::Text] && this->textLogger == nullptr) {
      this->textLogger = std::make_unique<LogTextStream<TextLogger>>();
    }
    this->enabledLoggers = enabledLoggers;
  }

  MasterLogging &operator<<(std::ostream &(*f)(std::ostream &)) {
    if (onMaster) {
      if (this->enabledLoggers[LoggingBackend::Text]) *this->textLogger << f;
    }
    return *this;
  };

  MasterLogging &operator<<(std::ios_base &(*f)(std::ios_base &)) {
    if (onMaster) {
      if (this->enabledLoggers[LoggingBackend::Text]) *this->textLogger << f;
    }

    return *this;
  };

  void flush() {
    if (!onMaster) return;
    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->flush();
    if (this->enabledLoggers[LoggingBackend::JSON]) this->jsonLogger->flush();
  }

  void AddErrorMsg(const std::string &msg) {
    if (!onMaster) return;
    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->AddErrorMsg(msg);
  }

  void AddWarningMsg(const std::string &msg) {
    if (!onMaster) return;
    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->AddWarningMsg(msg);
  }

  void CommunicateWarnings() {
    if (!onMaster) return;
    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->CommunicateWarnings();
  }

  void StartBlock(const std::string &blockName = "") {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->StartBlock(blockName);
    if (this->enabledLoggers[LoggingBackend::JSON]) this->jsonLogger->StartBlock(blockName);
  }

  void EndBlock(bool mute = false, const std::string &msg = "") {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->EndBlock(mute, msg);
    if (this->enabledLoggers[LoggingBackend::JSON]) this->jsonLogger->EndBlock();
  }

  void EndBlock(int classVerbose, int verbose = 0) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->EndBlock(classVerbose, verbose);
    if (this->enabledLoggers[LoggingBackend::JSON]) this->jsonLogger->EndBlock();
  }

  void StartBlockWithoutIndent(const std::string &blockName = "") {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->StartBlockWithoutIndent(blockName);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->StartBlock(blockName); // JSON does not have any indentation
  }

  void EndBlockWithoutIndent(bool mute = false) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->EndBlockWithoutIndent(mute);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->EndBlock(); // JSON does not have any indentation;
  }

  void startPoutBlock() {
    if (this->enabledLoggers[LoggingBackend::JSON]) { this->jsonLogger->startPoutBlock(); }
  }

  void endPoutBlock() {
    if (this->enabledLoggers[LoggingBackend::JSON]) { this->jsonLogger->EndBlock(); }
  }

  void printFromPLogging(const std::string &line_p, int proc) {
    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->printFromPLogging(line_p, proc);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->printFromPLogging(line_p, proc);
  }

  /// @brief Formats a key-value pair with padding
  /// @tparam T
  /// @param key
  /// @param val
  template<class T>
  void printEntryRaw(const std::string &key, const T &value) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text]) this->textLogger->printEntryRaw(key, value);
    if (this->enabledLoggers[LoggingBackend::JSON]) this->jsonLogger->printEntryRaw(key, value);
  }

  void PrintEntries(const std::string &title, int msgVerbose,
                    const std::map<std::string, std::vector<double>> &entries) {
    std::vector<PrintInfoEntry<std::vector<double>>> entriesVec;
    for (auto &[name, values] : entries) {
      entriesVec.push_back({name, values});
    }
    PrintInfo(title, msgVerbose, entriesVec);
  }

  template<typename... Types>
  void PrintEntries(const std::string &title, int msgVerbose, const Types &...types) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->PrintEntries(title, msgVerbose, types...);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->PrintEntries(title, msgVerbose, types...);
  }

  template<typename... Types>
  void PrintInfo(const std::string &className, int classVerbose, const Types &...types) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->PrintInfo(className, classVerbose, types...);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->PrintInfo(className, classVerbose, types...);
  }

  template<typename T>
  void PrintInfo(const std::string &className, int classVerbose,
                 const std::vector<PrintInfoEntry<T>> &entries) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->PrintInfo(className, classVerbose, entries);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->PrintInfo(className, classVerbose, entries);
  }

  // void PrintIteration(int classVerbose, const PrintIterEntry<double> &entries) {
  //   if (!onMaster)
  //     return;

  //   if (this->enabledLoggers[LoggingBackend::Text])
  //     this->textLogger->PrintIteration(classVerbose, entries);
  //   if (this->enabledLoggers[LoggingBackend::JSON])
  //     this->jsonLogger->PrintIteration(classVerbose, entries);
  // }

  template<typename... Types>
  void PrintIteration(int classVerbose, const PrintIterEntry<Types> &...types) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->PrintIteration(classVerbose, types...);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->PrintIteration(classVerbose, types...);
  }

  void PrintIteration(int classVerbose, const std::span<PrintIterEntry<double>> entries) {
    if (!onMaster) return;

    if (this->enabledLoggers[LoggingBackend::Text])
      this->textLogger->PrintIteration(classVerbose, entries);
    if (this->enabledLoggers[LoggingBackend::JSON])
      this->jsonLogger->PrintIteration(classVerbose, entries);
  }

  // All the types that can be handled by the backends should be handled by them
  // directly
  template<LogTextStreamable<TextLogger> T>
  MasterLogging &operator<<(const T value) {
    if (onMaster && this->enabledLoggers[LoggingBackend::Text]) *this->textLogger << value;

    return *this;
  };

  void PrintBuildInfo();

  // All processes tell the master about the environment they are running
  // on, so the master process can write that information to the logfile.
  //
  // DANGER: If any process fails to call this function on startup, master will
  //         hang.
  void CommunicateProcessEnvironment();
};

class ProcessLogging {
private:
  static LogTextStream<ProcessLogging> *instance;

  std::ostringstream buffer;
  bool doubleBlock = false;
  std::string doubleOutputFormat = "\%13.8f";
protected:
  ProcessLogging() {}
public:
  inline static LogTextStream<ProcessLogging> &the() {
    if (!instance) instance = new LogTextStream<ProcessLogging>();

    return *instance;
  }

  [[nodiscard]]
  int getPrecision() {
    return buffer.precision();
  }

  void setPrecision(const int precision) { buffer.precision(precision); }

  ProcessLogging &operator<<(std::ostream &(*f)(std::ostream &));

  ProcessLogging &operator<<(std::ios_base &(*f)(std::ios_base &)) {
    buffer << f;
    return *this;
  }

  ProcessLogging &operator<<(const double &);

  template<Ostreamable Streamable>
  ProcessLogging &operator<<(const Streamable streamable) {
    buffer << streamable;
    return *this;
  }
};

// Macros
#define mout MasterLogging::the()

#define vout(i)                                                                                    \
  if (verbose >= i) mout

#define tout(i)                                                                                    \
  if (TimeLevel > i) mout

#define OUT(s) __STRING(s) << " =" << std::endl << s

#define DOUT(s) __STRING(s) << ":" << s << " "

#define MOUT(s) mout << OUT(s) << std::endl;

#define pout ProcessLogging::the()

#define pvout(i)                                                                                   \
  if (verbose >= i) pout

#define POUT(s) pout << OUT(s) << std::endl;

#define vpout(i)                                                                                   \
  if (verbose >= i) pout

#define ppout(s)                                                                                   \
  if (PPM->proc() == s) pout

#endif