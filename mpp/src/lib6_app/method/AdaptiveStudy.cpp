#include "AdaptiveStudy.hpp"

std::shared_ptr<AdaptiveStudy> CreateAdaptiveStudy(PDESolverConfig conf, AdaptiveConfig aconf) {
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
      return std::make_unique<AdaptiveStudySTAcoustic>(conf, aconf);
    }
    if (conf.modelName == "STAcousticGL") {
      return std::make_unique<AdaptiveStudySTAcousticGL>(conf, aconf);
    }
    if (conf.modelName == "STAcousticMF") {
      return std::make_unique<AdaptiveStudySTAcousticMF>(conf, aconf);
    }
    if (conf.modelName == "STAcousticPG") {
      return std::make_unique<AdaptiveStudySTAcousticPG>(conf, aconf);
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
