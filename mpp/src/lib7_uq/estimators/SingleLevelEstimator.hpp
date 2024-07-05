#ifndef MLMC_MC_HPP
#define MLMC_MC_HPP

#include "MultiSampleFEM.hpp"
#include "WelfordAggregate.hpp"
#include "MeshesCreator.hpp"
#include <utility>


/*
 * A SingleLevelEstimator can represent a Monte Carlo,
 *    quasi-Monte Carlo or a Stochastic Collocation method
 *    depending on the sample generator used in IStochasticProblem
 */

struct SLEstimatorConfig {
  int initLevel = 0;

  int initSamples = 0;

  double epsilon = 0.0;

  bool parallel = true;

  bool onlyFine = false;

  bool delayCostUpdate = false;

  bool delayTotalUpdate = false;

  bool delayParallelUpdate = false;

  MultiSampleFEMConfig msFEMConfig;

  explicit SLEstimatorConfig() {
    if (!Config::IsInitialized()) return;

    Config::Get("Level", initLevel);
    Config::Get("epsilon", epsilon);
    Config::Get("OnlyFine", onlyFine);
    Config::Get("Samples", initSamples);
    Config::Get("ParallelEstimator", parallel);
    Config::Get("DelayCostUpdate", delayCostUpdate);
    Config::Get("DelayTotalUpdate", delayTotalUpdate);
    Config::Get("DelayParallelUpdate", delayParallelUpdate);
  }

  SLEstimatorConfig WithInitLevel(int level);

  SLEstimatorConfig WithInitSamples(int samples);

  SLEstimatorConfig WithEpsilon(double targetEpsilon);

  SLEstimatorConfig WithParallel(bool parallelEstimator);

  SLEstimatorConfig WithOnlyFine(bool onlyComputeOnFine);

  SLEstimatorConfig WithDelayTotalUpdate(bool delayTotUpdate);

  SLEstimatorConfig WithDelayParallelUpdate(bool delayParaUpdate);

  SLEstimatorConfig WithMultiSampleFEMConfig(const MultiSampleFEMConfig &conf);
};

std::ostream &operator<<(std::ostream &s, const SLEstimatorConfig &conf);

class SingleLevelEstimator {
private:
  int plotting = 0;

  int verbose = 1;

  int dM = 0;

  int dMcomm = 0;

  int commSplit = 0;

  WelfordAggregate aggregate;

  SLEstimatorConfig slEstmConf;

  std::unique_ptr<MultiSampleFEM> msFEM;

  const int &initLevel = slEstmConf.initLevel;

  void updateSamplesIfNeeded();

public:

  explicit SingleLevelEstimator(WelfordAggregate aggregate) :
      SingleLevelEstimator(SLEstimatorConfig(), std::move(aggregate)) {}

  explicit SingleLevelEstimator(const SLEstimatorConfig &slEstmConf) :
      SingleLevelEstimator(slEstmConf, WelfordAggregate()) {}

  SingleLevelEstimator(const SingleLevelEstimator &singleLevelEstimator) :
      SingleLevelEstimator(singleLevelEstimator.slEstmConf, singleLevelEstimator.aggregate) {}

  SingleLevelEstimator() : SingleLevelEstimator(SLEstimatorConfig(), WelfordAggregate()) {}

  SingleLevelEstimator(const SLEstimatorConfig &slEstmConf, WelfordAggregate aggregate) :
      slEstmConf(slEstmConf), aggregate(std::move(aggregate)) {

    UpdateSampleAmount(slEstmConf.initSamples);
    Config::Get("SLEstimatorVerbose", verbose);
    Config::Get("SLEstimatorPlotting", plotting);

    msFEM = CreateUniqueMSFEM(
        this->slEstmConf.msFEMConfig
            .WithCommSplit(commSplit)
            .WithCoarseLevel(slEstmConf.onlyFine ? initLevel : initLevel - 1)
            .WithFineLevel(initLevel)
    );
  }

  void Method();

  int CommSplit() const { return commSplit; }

  int SamplesToCompute() const { return dM; }

  double CellsOn(const LevelPair &meshIndex) const {
    return msFEM->GetMeshes()[meshIndex].CellCountGeometry();
//    * pow()
  }

  double Cost() const { return aggregate.cost; }

  double Value() const { return aggregate.Value(); }

  int Level() const { return slEstmConf.initLevel; }

  int SamplesToComputeOnComm() const { return dMcomm; }

  const MultiSampleFEM &GetMSFEM() const { return *msFEM; }

  double TotalError() const { return aggregate.errors.rmse; }

  double NumericError() const { return aggregate.errors.disc; }

  double StochasticError() const { return aggregate.errors.input; }

  const WelfordAggregate &GetAggregate() const { return aggregate; }

  void UpdateSVarYAndCostPerSample(double sVarY, double expectedCostPerSample);

  void UpdateOnComm(const SampleSolution &fSolution, const SampleSolution &cSolution);

  void UpdateSampleAmount(int requestedSamples);

  void UpdateCost(double newRoundCost);

  void UpdateParallel();

  void UpdateErrors();

  void UpdateTotal();

  int Index() const;

  void EstimatorResults() const {
    mout.PrintInfo("Estimator Results", verbose,
                   PrintInfoEntry("MSE", aggregate.errors.mse),
                   PrintInfoEntry("RMSE", aggregate.errors.rmse),
                   PrintInfoEntry("Cost", aggregate.cost),
                   PrintInfoEntry("Value", aggregate.Value()),
                   PrintInfoEntry("Epsilon", slEstmConf.epsilon),
                   PrintInfoEntry("Processes", PPM->Size(0)),
                   PrintInfoEntry("Disc Error", aggregate.errors.disc),
                   PrintInfoEntry("Input Error", aggregate.errors.input),
                   PrintInfoEntry("Used CPU sec", aggregate.cost * PPM->Size(0))
    );
  };

  static std::string Name() { return "SingleLevelEstimator"; }
};

#endif //MLMC_MC_HPP