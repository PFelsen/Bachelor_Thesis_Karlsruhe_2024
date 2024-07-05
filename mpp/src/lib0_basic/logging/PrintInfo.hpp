#include <string>
#include <vector>

#ifndef PRINTINFO_HPP
#define PRINTINFO_HPP

struct SeparatorClass {};

template<typename T>
class PrintInfoEntry {
  std::string msg;
  T value;
  int verbose;
public:
  PrintInfoEntry(const std::string &msg, const T &value, int verbose = 1) :
      msg(msg), value(value), verbose(verbose) {}

  const std::string &Message() const { return msg; }

  const T &Value() const { return value; }

  int Verbose() const { return verbose; }

  PrintInfoEntry &operator=(const PrintInfoEntry &entry) {
    msg = entry.msg;
    value = entry.value;
    verbose = entry.verbose;
    return *this;
  }

  [[nodiscard]]
  bool isPrintable(int maxVerbosityToDisplay) const noexcept {
    return verbose <= maxVerbosityToDisplay;
  }
};

/// @return true if any of the entries will be displayed at the specified
/// verbosity, false otherwise
template<typename... Types>
[[nodiscard]]
constexpr bool hasPrintableEntries(int maxVerbosityToDisplay,
                                   const PrintInfoEntry<Types>... entries) noexcept {
  return (... || entries.isPrintable(maxVerbosityToDisplay));
}

/// @return true if any of the entries will be displayed at the specified verbosity, false otherwise
template<typename T>
[[nodiscard]]
bool hasPrintableEntries(int maxVerbosityToDisplay,
                         const std::vector<PrintInfoEntry<T>> &entries) noexcept {
  for (const auto &entry : entries) {
    if (entry.isPrintable(maxVerbosityToDisplay)) { return true; }
  }
  return false;
}

template<typename T>
bool operator<(const PrintInfoEntry<T> &e1, const PrintInfoEntry<T> &e2) {
  std::string message1 = e1.Message();
  std::string message2 = e2.Message();
  std::transform(message1.begin(), message1.end(), message1.begin(),
                 [](auto character) { return std::tolower(character); });
  std::transform(message2.begin(), message2.end(), message2.begin(),
                 [](auto character) { return std::tolower(character); });

  return message1 < message2;
}

template<>
class PrintInfoEntry<SeparatorClass> {
  int verbose;
  int length;
public:
  PrintInfoEntry(int length = 50, int verbose = 1) : verbose(verbose), length(length) {}

  int Verbose() const { return verbose; }

  int Lenght() const { return length; }
};

template<typename T>
class PrintIterEntry {
  const std::string msg;
  const T value;
  int verbose;
  int setw;
public:
  PrintIterEntry(const std::string &msg, const T &value, int setw, int verbose = 1) :
      msg(msg), value(value), verbose(verbose), setw(setw) {}

  const std::string &Message() const { return msg; }

  const T &Value() const { return value; }

  int Verbose() const { return verbose; }

  int SetW() const { return setw; }

  bool isPrintable(int classVerbose) { return classVerbose >= this.Verbose(); }

  template<typename... Types>
  static bool hasPrintableEntries(int classVerbose, const PrintIterEntry<T> &entry,
                                  const Types &...types) {
    if (entry.isPrintable(classVerbose)) { return true; }
    if (sizeof...(types) == 0) { return false; }
    return hasPrintableEntries(classVerbose, types...);
  }
};

typedef PrintInfoEntry<std::string> InfoEntry;

typedef std::vector<InfoEntry> InfoEntries;

typedef PrintInfoEntry<SeparatorClass> PrintInfoSeparator;

#endif