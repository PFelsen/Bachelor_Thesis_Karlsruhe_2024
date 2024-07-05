#ifndef RESULTS_HPP
#define RESULTS_HPP


#include "Assertion.hpp"
#include "Logging.hpp"
#include "Functools.hpp"
#include "PrintUtil.hpp"
#include "StringUtil.hpp"

#include <string>
#include <map>
#include <variant>
#include <vector>



using ResultItem = std::variant<double, int, size_t, std::string>;

inline std::ostream &operator<<(std::ostream &l, ResultItem item) {
  auto log = ft::Overload{
      [&l, item](double) { l << std::get<double>(item); },
      [&l, item](int) { l << std::get<int>(item); },
      [&l, item](std::string &) { l << std::get<std::string>(item); },
      [&l, item](size_t &) { l << std::get<std::size_t>(item);}
  };
  std::visit(log, item);
  return l;
}

class ResultList {
  std::vector<ResultItem> data;

public:
  void append(const ResultItem &item) {
    data.push_back(item);
  }

  ResultItem &last() {
    if (data.size() == 0) {
      THROW("ResultList not filled.")
    }
    return data[data.size() - 1];
  }

  template<class T>
  T last() {
    if (data.size() == 0) {
      THROW("ResultList not filled.")
    }
    if (!std::holds_alternative<T>(data[data.size() - 1])) {
      THROW("Wrong in ResultList.")
    }
    return std::get<T>(data[data.size() - 1]);
  }

  template<class Type>
  std::vector<Type> asType() {
    std::vector<Type> returnValue;
    for (ResultItem &ri : data) {
      if (std::holds_alternative<Type>(ri)) {
        returnValue.push_back(std::get<Type>(ri));
      } else {
        Exit("Wrong type!")
      }
    }
    return returnValue;
  }

  std::vector<ResultItem> asRaw() {
    return data;
  }
};

class Results {
private:
  int verbose = 0;
  std::map<std::string, ResultList> results;

public:

  Results() {
    Config::Get("ResultsVerbose", verbose);
  }

  void append(const std::string &name, const ResultItem &value);

  ResultList &operator[](const std::string &name);

  template<typename Type>
  std::vector<Type> get(const std::string &name) {
    auto ls = results.find(name);
    if (ls == results.end()) {
      Exit("name " + name + " not found in results.")
    }
    return ls->second.asType<Type>();
  }

  void clear() {
    results.clear();
  }


  std::string ToJON() {
    std::stringstream ss;
    ss.precision(16);
    ss << std::fixed;

    ss << "{" << endl;
    for (auto &[key, values]: results) {
      ss << "\"";
      ss << key;
      ss << "\": [";
      for (int i = 0; i < values.asRaw().size(); i++) {
        ss << values.asRaw()[i];
        if (i < values.asRaw().size() - 1) {
          ss << ",";
        }
      }
      ss << "]" << endl;
    }
    ss << "}";
    return ss.str();
  }

  void PrintInfo() {
    if (!PPM->Master(0)) {
      return;
    }
    for (auto &[key, value]: results) {
      if (std::holds_alternative<std::string>(value.last())) {
        vout(1).printEntryRaw(key, concat(value.asType<std::string>()));
      }
    }

    for (auto &[key, value]: results) {
      if (std::holds_alternative<int>(value.last())) {
        vout(1).printEntryRaw(key, to_string(value.asType<int>()));
      }
      if (std::holds_alternative<std::size_t>(value.last())) {
        vout(1).printEntryRaw(key, to_string(value.asType<std::size_t>()));
      }
    }
    for (auto &[key, value]: results) {
      if (std::holds_alternative<double>(value.last())) {
        PrintValues(key, value.asType<double>(), 100, verbose);
      }
    }
  }
};

#endif //RESULTS_HPP
