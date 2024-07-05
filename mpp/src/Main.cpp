#include "Main.hpp"
#include "m++.hpp"


int main(int argv, char **argc) {
  Config::SetConfPath(std::string(ProjectMppDir) + "/conf/");
  Mpp::initialize(&argv, argc);

  Config::PrintInfo();

  std::string run = "single";
  Config::Get("run", run);

  if (run == "single") SingleExperiment();
  else if (run == "convergence") ConvergenceExperiment();
  else if (run == "adaptive") AdaptiveExperiment();
  else if (run == "adaptive_convergence") ; // TODO

#if USE_SPACETIME
  else if (run == "FWIST") FWIExperimentST();
#endif

#if BUILD_UQ
  else if (run == "SGD") SGDExperiment();
  else if (run == "MLMC") MLMCExperiment();
  else if (run == "MC") MCExperiment();
  else if (run == "SMC") SMCExperiment();
#endif
  return 0;
}
