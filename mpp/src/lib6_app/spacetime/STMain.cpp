#include "m++.hpp"
//#include "SpaceTime.hpp"
#include "STAdaptiveMethod.hpp"
#include "STSingleRunMethod.hpp"
#include "STConvergenceMethod.hpp"
#include "STAdaptiveConvergenceMethod.hpp"
#include "iostream"

int main(int argv, char **argc) {
  //std::cout << "STMain.cpp" << endl;
  //Config::PrintInfo();
  Config::SetConfPath(std::string(ProjectMppDir) + "/conf/");
  Mpp::initialize(&argv, argc);
  string run = "none";
  Config::Get("run", run);

  if (run == "single") {
    STSingleRunMethod().run();
  } else if (run == "adaptive") {
    STAdaptiveMethod().run();
  } else if (run == "convergence") {
    STConvergenceMethod().run();
  } else if (run == "adaptive_convergence") {
    STAdaptiveConvergenceMethod().run();
  }

  return 0;
}
