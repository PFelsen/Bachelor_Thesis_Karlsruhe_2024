
#include "Results.hpp"





void Results::append(const std::string &name, const ResultItem &value) {
  (*this)[name].append(value);
}

ResultList& Results::operator[](const std::string &name){
  auto it = results.find(name);
  if(it == results.end()){
    it = results.emplace(name, ResultList()).first;
  }
  return it->second;
}