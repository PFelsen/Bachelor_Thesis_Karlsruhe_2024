#ifndef LOGGING_BACKENDS_HPP
#define LOGGING_BACKENDS_HPP

#include <chrono>
#include <iostream>
#include <stack>
#include <string>

#include <json.hpp>

#include "M_IOFiles.hpp"
#include "PrintInfo.hpp"
#include "StringUtil.hpp"
#include "TimeDate.hpp"

#ifdef BUILD_IA
#include "IACInterval.hpp"
#include "IAInterval.hpp"
#endif

using json = nlohmann::json;

namespace std {
template<class CharT, class Traits>
std::basic_ostream<CharT, Traits> &beginl(std::basic_ostream<CharT, Traits> &os) {
  return os << std::flush << "\n";
}

template<class CharT, class Traits>
std::basic_ostream<CharT, Traits> &beginD(std::basic_ostream<CharT, Traits> &os) {
  return os;
}

template<class CharT, class Traits>
std::basic_ostream<CharT, Traits> &endD(std::basic_ostream<CharT, Traits> &os) {
  return os;
}
} // namespace std

using std::beginD;
using std::beginl;
using std::endD;
using std::endl;
using std::flush;

inline void to_json(json &j, const Date &d) {
  char *p = ctime(&d.t);
  p[strcspn(p, "\n")] = '\0';
  j = p;
}

inline void to_json(json &j, const Time &t) {
  int M = t.Minutes();
  double S = t.Seconds();
  std::string out = std::to_string(S);
  if (M == 0) {
    j = out + " seconds";
    return;
  }
  int H = t.Hours();
  S -= 60 * M;
  out = std::to_string(M) + ":" + std::to_string(S);
  if (H == 0) {
    j = out + " minutes";
    return;
  }
  M -= 60 * H;
  out = std::to_string(H) + ":" + std::to_string(M) + ":" + std::to_string(S);
  j = out + " hours";
}

/// @brief PrintInfoEntries that can be converted to json values
template<typename T>
concept JsonConvertible = requires(const T v) { json{v}; };

// Extends the API of an ostream-like type to include shift operators
// for container types like std::tuple and std::vector
//
// This allows for easy abstractions over TextLogger/Pout, which share
// a lot of functionality
template<typename S>
class LogTextStream : public S {
public:
  using S::operator<<;

  template<TupleLike TupleType>
    requires(!Ostreamable<TupleType>)
  LogTextStream<S> &operator<<(const TupleType &tuple) {
    std::apply([&](auto &...args) { ((*this << args), ...); }, tuple);

    return *this;
  };

  template<SizableIterableLike IterableType>
    requires(!Ostreamable<IterableType> && !MapIterableLike<IterableType>)
  LogTextStream<S> &operator<<(const IterableType &iterator) {
    for (const auto &element : iterator) {
      *this << element << std::endl;
    }
    return *this;
  };

  template<MapIterableLike MapType>
    requires(!Ostreamable<MapType>)
  LogTextStream<S> &operator<<(const MapType &map) {
    for (const auto &element : map) {
      *this << element.first << ": " << element.second << std::endl;
    }
    return *this;
  };
};

template<typename T, typename S>
concept LogTextStreamable =
    requires(LogTextStream<S> &logTextStream, T value) { logTextStream << value; };

struct TextBlock {
  std::chrono::time_point<std::chrono::system_clock> startTime;
  std::optional<std::string> name;
};

template<typename T>
std::string vec2str(std::vector<T> vec);

template<typename T>
std::string vec2str(std::vector<std::vector<T>> vec);

class TextLogger {
private:
  M_ofstream fileOut;

  bool fileEnabled = true;
  bool screenEnabled = true;

  static constexpr auto defaultIndent = "  ";
  static constexpr int defaultIndentLength = 2;
  Date startTime;
  std::stack<TextBlock> openBlocks;

  std::string logPrefix;
  std::string indent = "";

  std::string error_msg = "";
  int error_cnt = 0;
  std::string warning_msg = "";
  int warning_cnt = 0;

  bool lineEmpty = true;
  bool doubleBlock = false;
  bool onMaster = false;
  int numProcs;
  std::string doubleOutputFormat = "\%13.8f";

  std::stack<Date> timeStamps;
protected:
  void DecreaseIndent() { indent = indent.erase(0, defaultIndentLength); }

  void IncreaseIndent() { indent = indent.append(defaultIndent); }

  // Internal utility to manage writing to all enabled output streams
  // (console/file/both)
  template<typename T>
  inline TextLogger &write(T value) {
    if (fileEnabled) fileOut << value;
    if (screenEnabled) std::cout << value;
    return *this;
  }

  void printIndentMsg() {
    if (!lineEmpty) return;

    std::string msg = logPrefix + indent;
    if (!openBlocks.empty() && openBlocks.top().name.has_value())
      msg += openBlocks.top().name.value() + ": ";

    this->write(msg);

    lineEmpty = false;
  }

  [[nodiscard]]
  int getPrecision() {
    return fileOut.precision();
  }

  void setPrecision(int precision) {
    if (fileEnabled) fileOut.precision(precision);
    if (screenEnabled) std::cout.precision(precision);
  }
public:
  TextLogger();

  ~TextLogger();

  void setFileEnabled(bool to) { this->fileEnabled = to; }

  void setScreenEnabled(bool to) { this->screenEnabled = to; }

  TextLogger &operator<<(std::ostream &(*f)(std::ostream &));

  void StartBlock(const std::string &blockName = "");

  void StartPrintInfoBlock(const std::string &blockName);

  void EndBlock(bool mute = false, const std::string &msg = "");

  void EndBlock(int classVerbose, int verbose = 0);

  void EndPrintInfoBlock() {
    DecreaseIndent();
    this->endl();
  }

  void StartBlockWithoutIndent(const std::string &blockName = "");

  void EndBlockWithoutIndent(bool mute = false);

  void printFromPLogging(const std::string &line_p, int proc);

  void AddErrorMsg(const std::string &msg);

  void AddWarningMsg(const std::string &msg);

  void CommunicateWarnings(){};

  void flush() {
    if (fileEnabled) fileOut << std::flush;
    if (screenEnabled) std::cout << std::flush;
  }

  void endl() { *this << std::endl; }

  void beginl() {
    if (fileEnabled) fileOut << std::beginl;
    if (screenEnabled) std::cout << std::beginl;
  }

  TextLogger &operator<<(std::ios_base &(*f)(std::ios_base &)) {
    this->write(f);
    return *this;
  }

  TextLogger &operator<<(std::_Setw f) {
    if (fileEnabled) fileOut.width(f._M_n);
    if (screenEnabled) std::cout.width(f._M_n);
    return *this;
  }

  TextLogger &operator<<(const std::chrono::duration<double> duration) {
    printIndentMsg();
    this->write(std::to_string(duration.count()) + " seconds");
    return *this;
  }

  TextLogger &operator<<(const std::chrono::time_point<std::chrono::system_clock> &point) {
    printIndentMsg();
    std::stringstream format;
    const std::time_t time = std::chrono::system_clock::to_time_t(point);
    std::tm localTime = *std::localtime(&time);
    format << std::put_time(&localTime, "%a %b %d %H:%M:%S %Y");
    const auto formatedString = format.str();
    this->write(formatedString.c_str());
    return *this;
  }

  TextLogger &operator<<(const std::string &s) {
    std::vector<std::string> parts = splitWithoutTrim(s, "\n");
    printIndentMsg();
    this->write(parts[0]);
    for (const auto &msg : parts | std::views::drop(1)) {
      *this << std::endl;
      printIndentMsg();
      this->write(msg);
    }
    return *this;
  }

  TextLogger &operator<<(const char *s) {
    std::string msg(s);
    return *this << msg;
  }

  TextLogger &operator<<(const double &s);

  template<class T>
  inline void printEntryRaw(const std::string &key, const T &val) {
    std::string msg(key + ": ");
    int l = 40 - int(key.size());
    if (l > 0) msg.append(std::string(l, '.')).append(" ");
    (*this) << msg;
    if constexpr (std::is_same<T, std::vector<int>>::value
                  || std::is_same<T, std::vector<double>>::value
                  || std::is_same<T, std::vector<std::string>>::value
                  || std::is_same<T, std::vector<std::vector<int>>>::value
                  || std::is_same<T, std::vector<std::vector<double>>>::value
                  || std::is_same<T, std::vector<std::vector<std::string>>>::value) {
      *(LogTextStream<TextLogger> *)(this) << vec2str(val);
    } else {
      *(LogTextStream<TextLogger> *)(this) << val;
    }
    this->endl();
  }

  template<Ostreamable Streamable>
  TextLogger &operator<<(const Streamable streamable) {
    printIndentMsg();
    if (screenEnabled) std::cout << streamable;
    if (fileEnabled) fileOut << streamable;
    return *this;
  }

  template<class T>
  void printIterEntryRaw(const std::string &msg, const T &val, int setw) {
    // TODO vereinheitlichen mit printEntry und width übergeben als
    // optional parameter
    (*this) << msg << "=" << std::setw(setw) << std::left << val << " ";
  }

  template<typename T>
  void printOldEntry(int maxVerbosityToDisplay, const PrintInfoEntry<T> &entry) {
    if (maxVerbosityToDisplay < entry.Verbose()) return;
    printEntryRaw(entry.Message(), entry.Value());
  }

  template<typename... Types>
  void PrintEntries(const std::string &title, int maxVerbosityToDisplay,
                    const PrintInfoEntry<Types>... entries) {
    if (!hasPrintableEntries(maxVerbosityToDisplay, entries...)) return;

    StartPrintInfoBlock(title);
    (printOldEntry(maxVerbosityToDisplay, entries), ...);
    EndPrintInfoBlock();
  }

  template<typename... Types>
  void PrintInfo(const std::string &title, int maxVerbosityToDisplay, const Types &...types) {
    PrintEntries(title + " Info", maxVerbosityToDisplay, types...);
  }

  template<typename T>
  void PrintInfo(const std::string &title, int maxVerbosityToDisplay,
                 const std::vector<PrintInfoEntry<T>> &entries) {
    if (!hasPrintableEntries(maxVerbosityToDisplay, entries)) return;

    StartPrintInfoBlock(title + " Info");
    for (const auto &entry : entries)
      printOldEntry(maxVerbosityToDisplay, entry);
    EndPrintInfoBlock();
  };

  template<typename T>
  void PrintIteration(int maxVerbosityToDisplay, const std::span<PrintIterEntry<T>> entries) {
    for (const auto &entry : entries) {
      printIterEntryRaw(entry.Message(), entry.Value(), entry.SetW());
    }
    *this << std::endl;
  }

  template<typename... Types>
  void PrintIteration(int maxVerbosityToDisplay, const Types &...types) {
    (printIterEntryRaw(types.Message(), types.Value(), types.SetW()), ...);
    *this << std::endl;
  }
};

struct JsonBlock {
  // The blocks data as a json object
  json data;
  // The JSON Logger ignores blocks without a name.
  // To track the corresponding endBlock() calls, it is necessary
  // to count how many of these "unnamed blocks" were opened since
  // this block was started, as that many endBlock() calls need to be ignored
  int unnamedOpenGroups;

  // How many pout blocks were opened in this block
  int numPoutBlocks;

  std::string blockName;

  std::chrono::time_point<std::chrono::system_clock> startTime;

  /// @brief Whether or not closing this block will terminate a printIterationContext
  bool terminatesPrintIteration = false;
};

class JsonLogger {
private:
  bool dirty;
  std::deque<JsonBlock> openBlocks;
  M_ofstream fileOut;
  bool onMaster;
  std::optional<json::array_t> printIterationContext;
public:
  JsonLogger();

  ~JsonLogger();

  void StartBlock(const std::string &blockName);

  void EndBlock();

  json::array_t &startPrintIterationContext() {
    // If we are not inside a printIterationContext, start one
    if (!printIterationContext.has_value()) {
      printIterationContext = json::array();
      openBlocks.back().terminatesPrintIteration = true;
    }

    return printIterationContext.value();
  }

  void endPrintIterationContext() {
    // If we are not inside a printIterationContext then there's nothing to do
    if (!printIterationContext.has_value() || printIterationContext.value().empty()) return;
    auto context = printIterationContext.value();

    const auto start = std::chrono::system_clock::now();
    std::string blockName = "iteration_" + std::to_string(this->openBlocks.back().numPoutBlocks);
    openBlocks.back().numPoutBlocks++;

    openBlocks.back().data[blockName] = json(context);
    printIterationContext.reset();
  }

  void startPoutBlock() {
    const auto start = std::chrono::system_clock::now();
    std::string blockName = "pout_" + std::to_string(this->openBlocks.back().numPoutBlocks);
    this->openBlocks.back().numPoutBlocks++;

    json block(json::value_t::array);
    openBlocks.push_back(JsonBlock{block, 0, 0, blockName, start});
  }

  void printFromPLogging(const std::string &line_p, int proc) {
    openBlocks.back().data.push_back(json(line_p));
    dirty = true;
  }

  void flush();

  template<typename T>
  void printEntryRaw(const std::string &key, const T &value) {
    if constexpr (requires { json(value); }) {
      openBlocks.back().data[key] = json(value);
      dirty = true;
    }
  }

  template<typename T>
  void printOldEntry(int maxVerbosityToDisplay, const PrintInfoEntry<T> &entry) {
    if (maxVerbosityToDisplay < entry.Verbose()) return;
    printEntryRaw(entry.Message(), entry.Value());
  }

  template<typename... Types>
  void PrintEntries(const std::string &title, int maxVerbosityToDisplay,
                    const PrintInfoEntry<Types> &...entries) {
    if (!hasPrintableEntries(maxVerbosityToDisplay, entries...)) return;

    endPrintIterationContext();
    StartBlock(title);
    (printOldEntry(maxVerbosityToDisplay, entries), ...);
    EndBlock();
  }

  template<typename... Types>
  void PrintInfo(const std::string &title, int maxVerbosityToDisplay, const Types &...types) {
    PrintEntries(title + " Info", maxVerbosityToDisplay, types...);
  }

  template<typename T>
  void PrintInfo(const std::string &title, int maxVerbosityToDisplay,
                 const std::vector<PrintInfoEntry<T>> &entries) {
    StartBlock(title + " Info");
    for (const auto &entry : entries)
      printOldEntry(maxVerbosityToDisplay, entry);
    EndBlock();
  };

  template<typename T>
  void PrintIteration(int maxVerbosityToDisplay, const std::span<PrintIterEntry<T>> entries) {
    startPrintIterationContext();

    json iteration_value = {};
    for (const auto &entry : entries) {
      iteration_value[entry.Message()] = entry.Value();
    }
    printIterationContext->push_back(iteration_value);
  }

  template<JsonConvertible... Types>
  void PrintIteration(int msgVerbose, const PrintIterEntry<Types> &...types) {
    json iteration_value = {{types.Message(), types.Value()}...};

    startPrintIterationContext();
    printIterationContext->push_back(iteration_value);
  }
};

#endif