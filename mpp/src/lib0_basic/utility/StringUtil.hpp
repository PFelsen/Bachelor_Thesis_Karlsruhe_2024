#ifndef _STRINGUTIL_CPP_
#define _STRINGUTIL_CPP_

#include <sstream>
#include <string>
#include <vector>

// trim from start (in place)
void ltrim(std::string &s);

// trim from end (in place)
static inline void rtrim(std::string &s);

// trim from both ends (in place)
void trim(std::string &s);

// trim from start (copying)
std::string ltrim_copy(std::string s);

// trim from end (copying)
std::string rtrim_copy(std::string s);

// trim from both ends (copying)
std::string trim_copy(std::string s);

std::vector<std::string> split(std::string toSplit, std::string delimiters);

std::vector<std::string> splitWithoutTrim(std::string toSplit, std::string delimiters);

std::vector<int> splitToInt(std::string toSplit, std::string delimiters);

std::vector<double> splitToDouble(std::string toSplit, std::string delimiters);

bool contains(const std::string &s, const std::string &contained);

std::string join(const std::string &delimiter, const std::vector<std::string> &strs);

std::string lowerCase(const std::string &str);

std::string upperCase(const std::string &str);

std::string capitalize(const std::string &str);

std::string makeStringPathLike(const std::string &str);

std::string removeBrackets(const std::string &value);

template<typename T>
std::string concat(std::vector<T> strs) {
  std::stringstream ss;
  ss << "[ ";
  if (strs.empty()) {
    ss << "]";
    return ss.str();
  }
  ss << strs[0];
  for (int i = 1; i < strs.size(); i++) {
    ss << ", " << strs[i];
  }
  ss << " ]";
  return ss.str();
}

template<typename T>
T parse(const std::string &value);

#endif // of #ifndef _STRINGUTIL_CPP_
