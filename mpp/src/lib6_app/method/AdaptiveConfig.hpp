#ifndef MPP_ADAPTIVECONFIG_H
#define MPP_ADAPTIVECONFIG_H

#include "Config.hpp"
#include <string>

struct AdaptiveData {
public:
  int refinement_steps = -1;
  std::string errorEstimator = "none";
  double theta = -1;
  double theta_factor = 1; // 1 is neutral factor
  double theta_min = -1;
  std::string refine_by = "none";


  void ReadConfig() {
    if (errorEstimator == "none")
      Config::Get("ErrorEstimator", errorEstimator);
    if (refinement_steps == -1) {
      Config::Get("refinement_steps", refinement_steps);
    }
    Config::Get("theta", theta);
    Config::Get("theta_factor", theta_factor);
    Config::Get("theta_min", theta_min);
    Config::Get("refine_by", refine_by);
  }

};

struct AdaptiveConfig {
public:
  AdaptiveData data;

  AdaptiveConfig(AdaptiveData data) : data(data) {
    data.ReadConfig();
  }

  AdaptiveConfig() {
    data.ReadConfig();
  }
};


#endif //MPP_ADAPTIVECONFIG_H
