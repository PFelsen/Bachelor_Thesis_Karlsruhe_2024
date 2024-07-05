#ifndef MPP_ADAPTIVESTUDY_HPP
#define MPP_ADAPTIVESTUDY_HPP

#include "ConvergenceStudy.hpp"
#include "AdaptiveConfig.hpp"

class AdaptiveStudy : public ConvergenceStudy {
private:
  vector<LevelPair> GetLevelsForAdaptivity(LevelPair level) const {
    std::vector<LevelPair> levels{level};
    LevelPair current = level;
    do {
      current = current.NextInAdaptivity();
      levels.push_back(current);
    } while (current.adaptivityLevel != aconf.data.refinement_steps);
    return levels;
  }

  vector<LevelPair> getLevels() override {
    return GetLevelsForAdaptivity(getMeshes().FineLevel());
  }


public:
  const AdaptiveConfig aconf;

  explicit AdaptiveStudy(PDESolverConfig conf, AdaptiveConfig aconf)
  : ConvergenceStudy(std::move(conf)), aconf(std::move(aconf)) {
  }
};


#if USE_SPACETIME

#include "STAcousticMain.hpp"

template<class SpaceTimeAcousticPDESolver>
class AdaptiveStudySTAcousticT : public AdaptiveStudy {
  SpaceTimeAcousticPDESolver pdeSolver;
  std::shared_ptr<AcousticProblem> problem;

  const Meshes &getMeshes() override {
    return problem->GetMeshes();
  }

  Solution run(LevelPair level) override {
    Solution solution = pdeSolver.Run(problem, level);
    pdeSolver.EstimateErrorAndApplyAdaptivity(solution.vector, aconf);
    return solution;
  }

public:
  AdaptiveStudySTAcousticT(PDESolverConfig conf, AdaptiveConfig aconf) : AdaptiveStudy(conf, aconf), pdeSolver(conf) {
    problem = CreateAcousticProblemShared(conf.problemName);
    problem->CreateMeshes(MeshesCreator(conf.mesh).WithLevel(conf.level).WithPLevel(conf.pLevel));
  }
};

using AdaptiveStudySTAcoustic = AdaptiveStudySTAcousticT<STAcousticMain>;
using AdaptiveStudySTAcousticPG = AdaptiveStudySTAcousticT<STAcousticMainPG>;
using AdaptiveStudySTAcousticMF = AdaptiveStudySTAcousticT<STAcousticMainMF>;
using AdaptiveStudySTAcousticGL = AdaptiveStudySTAcousticT<STAcousticMainGL>;
#endif

std::shared_ptr<AdaptiveStudy> CreateAdaptiveStudy(PDESolverConfig conf, AdaptiveConfig aconf);

#endif
