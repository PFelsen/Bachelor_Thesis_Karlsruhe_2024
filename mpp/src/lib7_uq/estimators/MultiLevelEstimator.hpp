#ifndef MULTILEVELESTIMATOR_HPP
#define MULTILEVELESTIMATOR_HPP

#include "EstimatorMap.hpp"

struct ContinuationData {
  std::vector<double> costs{};

  std::vector<double> memory{};

  std::vector<double> values{};

  std::vector<double> epsilons{};

  std::vector<double> msErrors{};

  std::vector<double> rmsErrors{};

  std::vector<double> discErrors{};

  std::vector<double> disc2Errors{};

  std::vector<double> inputErrors{};

  std::vector<std::vector<int>> commSplits{};

  void PushBack(const EstimatorMap &estMap, double epsilon) {
    memory.push_back(MemoryLogger::TotalMemoryUsage());
    epsilons.push_back(epsilon);
    commSplits.push_back(estMap.GetCommSplitVector());
    inputErrors.push_back(estMap.GetErrors().input);
    disc2Errors.push_back(estMap.GetErrors().disc2);
    discErrors.push_back(estMap.GetErrors().disc);
    rmsErrors.push_back(estMap.GetErrors().rmse);
    msErrors.push_back(estMap.GetErrors().mse);
    values.push_back(estMap.Value());
    costs.push_back(estMap.Cost());
  }

  bool IsEmpty() const { return costs.empty(); }

  void PrintInfo() const;
};

class MultiLevelEstimator {
private:
  int plotting = 1;

  int verbose = 1;

  int attemptCtr = 0;

  double epsilon = 0.0;

  EstimatorMap estMap;

  ContinuationData contData;

  MLEstimatorConfig mlEstmConf;

  bool levelLoop();

  bool expectedToFinish();

  bool tryUpdateEpsilon();

  bool appendLevelIfNeeded();

  bool expectedToFitOnMemory();

  bool updateSamplesIfNeeded();

  bool timeBudgetIsExhausted();

  void updateData(std::map<int, std::pair<double, double>> timeMap);

public:
  explicit MultiLevelEstimator(MLEstimatorConfig mlEstmConf)
      : epsilon(mlEstmConf.epsilon), mlEstmConf(mlEstmConf),
        estMap(std::move(mlEstmConf)) {
    Config::Get("MLEstimatorPlotting", plotting);
    Config::Get("MLEstimatorVerbose", verbose);
  }

  ~MultiLevelEstimator() { estMap.clear(); }

  void MLE();

  void AMLE();

  void BMLE();

  void Method();

  void MultilevelResults() const;

  void ContinuationResults() const;

  void ExponentResults() const;

  double Cost() const { return estMap.Cost(); }

  double Value() const { return estMap.Value(); }

  const EstimatorMap &GetEstimatorMap() const { return estMap; }

  static std::string Name() { return "MultiLevelEstimator"; }

  void EstimatorResults() const {
    mout.PrintInfo("Estimator Results", verbose,
                   PrintInfoEntry("MSE", estMap.GetErrors().mse),
                   PrintInfoEntry("RMSE", estMap.GetErrors().rmse),
                   PrintInfoEntry("Cost", estMap.Cost()),
                   PrintInfoEntry("Value", estMap.Value()),
                   PrintInfoEntry("Epsilon", epsilon),
                   PrintInfoEntry("Processes", PPM->Size(0)),
                   PrintInfoEntry("Time Budget", mlEstmConf.timeBudget),
                   PrintInfoEntry("Cost Budget", mlEstmConf.timeBudget * PPM->Size(0)),
                   PrintInfoEntry("Disc Error", estMap.GetErrors().disc),
                   PrintInfoEntry("Disc2 Error", estMap.GetErrors().disc2),
                   PrintInfoEntry("Input Error", estMap.GetErrors().input),
                   PrintInfoEntry("Used CPU sec", estMap.Cost() * PPM->Size(0)),
                   PrintInfoEntry("Synchronization Efficiency",
                                  estMap.SynchronizationEfficiency()));
  };
};

#endif  // MULTILEVELESTIMATOR_HPP
