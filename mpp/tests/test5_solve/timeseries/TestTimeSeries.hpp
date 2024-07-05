#ifndef TESTTIMESERIES_HPP
#define TESTTIMESERIES_HPP

#include <map>
#include <string>


using ConfigMap = std::map<std::string, std::string>;

const ConfigMap testConfigMap = {{"startTime", "1.0"},
                                 {"endTime", "3.0"},
                                 {"stepSize", "2.0"},
                                 {"steps", "1"}};

#endif // TESTTIMESERIES_HPP
