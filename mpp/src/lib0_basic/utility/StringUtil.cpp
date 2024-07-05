
#include "StringUtil.hpp"
#include "Point.hpp"

#include <algorithm>
#include <cctype>
#include <locale>
#include <sstream>
#include <vector>

void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) { return !std::isspace(ch); }));
}

void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); }).base(),
          s.end());
}

std::string rtrim_copy(std::string s) {
  rtrim(s);
  return s;
}

void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

std::string trim_copy(std::string s) {
  trim(s);
  return s;
}

std::string ltrim_copy(std::string s) {
  ltrim(s);
  return s;
}

std::string join(const std::string &delimiter, const std::vector<std::string> &strs) {
  std::stringstream ss;
  if (strs.empty()) { return ss.str(); }
  ss << strs[0];
  for (int i = 1; i < strs.size(); i++) {
    ss << delimiter << strs[i];
  }
  return ss.str();
}

bool contains(const std::string &s, const std::string &contained) {
  return s.find(contained) != std::string::npos;
}

std::vector<std::string> splitWithoutTrim(std::string toSplit, std::string delimiters) {
  std::vector<std::string> out;
  size_t last = 0;
  size_t next = 0;
  toSplit += delimiters[0];
  while ((next = toSplit.find_first_of(delimiters, last)) != std::string::npos) {
    std::string line = toSplit.substr(last, next - last);
    last = next + 1;
    out.push_back(line);
  }
  return out;
}

std::vector<std::string> split(std::string toSplit, std::string delimiters) {
  std::vector<std::string> out;
  size_t last = 0;
  size_t next = 0;
  toSplit += delimiters[0];
  while ((next = toSplit.find_first_of(delimiters, last)) != std::string::npos) {
    std::string line = toSplit.substr(last, next - last);
    last = next + 1;
    trim(line);
    out.push_back(line);
  }
  return out;
}

std::vector<int> splitToInt(std::string toSplit, std::string delimiters) {
  trim(toSplit);
  auto splitted = split(toSplit, delimiters);

  std::vector<int> intSplitted;
  for (const auto &s : splitted) {
    if (!s.empty()) { intSplitted.emplace_back(std::stoi(s)); }
  }

  return intSplitted;
}

std::vector<double> splitToDouble(std::string toSplit, std::string delimiters) {
  trim(toSplit);
  auto splitted = split(toSplit, delimiters);

  std::vector<double> doubleSplitted;
  for (const auto &s : splitted) {
    if (!s.empty()) { doubleSplitted.emplace_back(std::stod(s)); }
  }

  return doubleSplitted;
}

std::string lowerCase(const std::string &str) {
  std::string lower(str);
  transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  return lower;
}

std::string upperCase(const std::string &str) {
  std::string upper(str);
  transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
  return upper;
}

std::string capitalize(const std::string &str) {
  std::string capStr(str);
  if (!str.empty()) {
    capStr[0] = std::toupper(str[0]);

    for (std::size_t i = 1; i < str.length(); ++i)
      capStr[i] = std::tolower(str[i]);
  }
  return capStr;
}

std::string makeStringPathLike(const std::string &str) {
  std::string path(str);
  if (path.back() != '/') { path.append("/"); }
  return path;
}

std::string removeBrackets(const std::string &value) {
  char brackets[] = "[](){}";
  std::string new_value = value;
  for (auto &bracket : brackets)
    new_value.erase(remove(new_value.begin(), new_value.end(), bracket), new_value.end());
  return new_value;
}

template<>
Point parse<Point>(const std::string &value) {
  std::vector<std::string> elems = split(value, ",");
  double coords[] = {0.0, 0.0, 0.0, 0.0};
  for (int index = 0; index < elems.size() && index < 4; index++) {
    std::string element = elems[index];
    if (!element.empty()) { coords[index] = stod(element); }
  }
  return Point(coords);
}

template<>
std::string parse<std::string>(const std::string &value) {
  return value;
}

template<>
bool parse<bool>(const std::string &value) {
  return (value == "true" || value == "True" || value != "0")
         && !(value == "false" || value == "False");
}

template<>
int parse<int>(const std::string &value) {
  return std::stoi(value);
}

template<>
size_t parse<size_t>(const std::string &value) {
  std::stringstream sstream(value);
  size_t result;
  sstream >> result;
  return result;
}

template<>
short parse<short>(const std::string &value) {
  return (short)std::stoi(value);
}

template<>
double parse<double>(const std::string &value) {
  return std::stod(value);
}

template<>
std::vector<int> parse<std::vector<int>>(const std::string &value) {
  std::string new_value = removeBrackets(value);
  std::vector<int> out;
  for (const std::string &element : split(new_value, ",")) {
    if (!element.empty()) { out.push_back(stoi(element)); }
  }
  return out;
}

template<>
std::vector<double> parse<std::vector<double>>(const std::string &value) {
  std::string new_value = removeBrackets(value);
  std::vector<double> out;
  for (const std::string &element : split(new_value, ",")) {
    if (!element.empty()) { out.push_back(stod(element)); }
  }
  return out;
}

template<>
std::vector<std::string> parse<std::vector<std::string>>(const std::string &value) {
  std::string new_value = removeBrackets(value);
  std::vector<std::string> out;
  for (const std::string &element : split(new_value, ",")) {
    if (!element.empty()) { out.push_back(element); }
  }
  return out;
}