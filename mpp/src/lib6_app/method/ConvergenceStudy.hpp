#ifndef MPP_CONVERGENCESTUDY_HPP
#define MPP_CONVERGENCESTUDY_HPP

#include "PDESolver.hpp"
#include "Matrix.hpp"
#include "Transfers.hpp"
#include "Plotting.hpp"


std::map<std::string, std::vector<double>>
convertResults(std::vector<std::map<std::string, double>> results);

class ConvergenceStudy {
private:
  PDESolverConfig conf;


  virtual vector<LevelPair> getLevels() {
    return GetLevels(getMeshes().PLevel(), getMeshes().FineLevel());
  }

  vector<LevelPair> GetLevels(LevelPair start, LevelPair end) {
    std::vector<LevelPair> levels{start};
    LevelPair current = start;
    do {
      current = LevelPair::Next(current, end);
      levels.push_back(current);
    } while (current != end);
    return levels;
  }

protected:
  std::vector<std::map<std::string, double>> results;

  virtual Solution run(LevelPair level) = 0;

  virtual const Meshes &getMeshes() = 0;

public:
  int verbose = 0;

  explicit ConvergenceStudy(PDESolverConfig conf) : conf(std::move(conf)) {
    Config::Get("ResultsVerbose", verbose);
  }

  std::vector<std::map<std::string, double>> GetResults() {
    return results;
  }

  virtual void Method() {
    for (LevelPair levels: getLevels()) {
      mout << "Convergence on Level " << levels.str() << endl;
      Solution solution = run(levels);
      results.push_back(solution.values);
      Plotting::Instance().Clear();
    }
    mout.PrintEntries("Convergence Results", verbose, convertResults(results));
  }
};

#if USE_SPACETIME

#include "STAcousticMain.hpp"

template<class SpaceTimeAcousticPDESolver>
class ConvergenceStudySTAcousticT : public ConvergenceStudy {
  SpaceTimeAcousticPDESolver pdeSolver;
  std::shared_ptr<AcousticProblem> problem;

  const Meshes &getMeshes() override {
    return problem->GetMeshes();
  }

  Solution run(LevelPair level) override {
    return pdeSolver.Run(problem, level);
  }

public:
  ConvergenceStudySTAcousticT(PDESolverConfig conf) : ConvergenceStudy(conf), pdeSolver(conf) {
    problem = CreateAcousticProblemShared(conf.problemName);
    problem->CreateMeshes(MeshesCreator(conf.mesh).WithLevel(conf.level).WithPLevel(conf.pLevel));
  }
};

using ConvergenceStudySTAcoustic = ConvergenceStudySTAcousticT<STAcousticMain>;
using ConvergenceStudySTAcousticPG = ConvergenceStudySTAcousticT<STAcousticMainPG>;
using ConvergenceStudySTAcousticMF = ConvergenceStudySTAcousticT<STAcousticMainMF>;
using ConvergenceStudySTAcousticGL = ConvergenceStudySTAcousticT<STAcousticMainGL>;
#endif

std::shared_ptr<ConvergenceStudy> CreateConvergenceStudy(PDESolverConfig conf);

#endif
