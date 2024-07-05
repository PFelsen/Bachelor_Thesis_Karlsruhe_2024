#include "ConvergenceStudy.hpp"

#if USE_SPACETIME
#include "STAcousticMain.hpp"
#endif


std::map<std::string, std::vector<double>>
convertResults(std::vector<std::map<std::string, double>> results) {
  std::map<std::string, std::vector<double>> res;
  for (auto &resultOnLevel: results) {
    for (auto &[name, value]: resultOnLevel) {
      res[name].push_back(value);
    }
  }
  return res;
}



std::shared_ptr<ConvergenceStudy> CreateConvergenceStudy(PDESolverConfig conf) {
  if (conf.modelName.find("Elliptic") != std::string::npos) {
    THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
  }

  if (conf.modelName.find("DGTransport") != std::string::npos) {
    THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
  }

  if (conf.modelName.find("Reaction") != std::string::npos) {
    THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
  }
#if USE_SPACETIME
  if (conf.modelName.find("STAcoustic") != std::string::npos) {
    if (conf.modelName == "STAcousticPDESolver") {
      THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
    }
    if (conf.modelName == "STAcoustic") {
      return std::make_unique<ConvergenceStudySTAcoustic>(conf);
    }
    if (conf.modelName == "STAcousticGL") {
      return std::make_unique<ConvergenceStudySTAcousticGL>(conf);
    }
    if (conf.modelName == "STAcousticMF") {
      return std::make_unique<ConvergenceStudySTAcousticMF>(conf);
    }
    if (conf.modelName == "STAcousticPG") {
      return std::make_unique<ConvergenceStudySTAcousticPG>(conf);
    }
  }
  if (conf.modelName == "STTransport") {
    THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
  }
#endif
  if (conf.modelName.find("VectorValued") != std::string::npos) {
    THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
  }

  if (conf.modelName.find("Acoustic") != std::string::npos) {
    THROW("implement CreateConvergenceStudy(PDESolverConfig conf) ")
  }

  Exit(conf.modelName + " not found")

}
