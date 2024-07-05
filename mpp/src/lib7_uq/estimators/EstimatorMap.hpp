#ifndef LEVELMAPS_HPP
#define LEVELMAPS_HPP

#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "SingleLevelEstimator.hpp"

struct MLEstimatorConfig {
  double eta = 0.9;

  double theta = 0.5;

  double epsilon = 0.0;

  double timeBudget = 0.0;

  double memoryBudget = 0.0;

  std::vector<int> initLevels{};

  std::vector<int> initSamples{};

  SLEstimatorConfig slEstmConf;

  explicit MLEstimatorConfig() {
    if (!Config::IsInitialized()) return;

    std::string wTimeStr = "00:00:00";
    Config::Get("initLevels", initLevels);
    Config::Get("initSamples", initSamples);
    Config::Get("MemoryBudget", memoryBudget);
    Config::Get("epsilon", epsilon);
    Config::Get("WTime", wTimeStr);
    Config::Get("theta", theta);
    Config::Get("eta", eta);

    double ratio_budget_reserved = 0.5;
    TRY {
      timeBudget = ratio_budget_reserved * (
          stoi(wTimeStr.substr(0, 2)) * 3600 +
          stoi(wTimeStr.substr(3, 2)) * 60 +
          stoi(wTimeStr.substr(6, 2)));
    }
    CATCH("Wrong time format")
  }

  MLEstimatorConfig WithSLEstimatorConfig(const SLEstimatorConfig &conf) {
    slEstmConf = conf;
    return *this;
  }

  MLEstimatorConfig WithInitSamples(const std::vector<int> &samples) {
    initSamples = samples;
    return *this;
  }

  MLEstimatorConfig WithInitLevel(const std::vector<int> &levels) {
    initLevels = levels;
    return *this;
  }

  MLEstimatorConfig WithTimeBudget(double timeBudgetAsDouble) {
    timeBudget = timeBudgetAsDouble;
    return *this;
  }

  MLEstimatorConfig WithTheta(double biasVarianceTradeOff) {
    theta = biasVarianceTradeOff;
    return *this;
  }

  MLEstimatorConfig WithEpsilon(double targetEpsilon) {
    epsilon = targetEpsilon;
    return *this;
  }

  MLEstimatorConfig WithEta(double errorReduction) {
    eta = errorReduction;
    return *this;
  }
};

std::ostream &operator<<(std::ostream &s, const MLEstimatorConfig &conf);


class EstimatorMap;

struct Exponents {
  double alpha = 0.0;

  double beta = 0.0;

  double gamma = 0.0;

  void ComputeAlpha(const EstimatorMap &mcMap);

  void ComputeBeta(const EstimatorMap &mcMap);

  void ComputeGamma(const EstimatorMap &mcMap);

  friend std::ostream &operator<<(std::ostream &s, const Exponents &exponents) {
    return s << "alpha=" << exponents.alpha << " beta=" << exponents.beta
             << " gamma=" << exponents.gamma;
  }
};

struct MultilevelErrors {
  double mse = 0.0;

  double rmse = 0.0;

  double disc = 0.0;

  double disc2 = 0.0;

  double input = 0.0;

  void EstimateErrors(const EstimatorMap &mcMap);

  static double EstimateNumeric(const EstimatorMap &mcMap);

  static double EstimateStochastic(const EstimatorMap &mcMap);

  friend std::ostream &operator<<(std::ostream &s,
                                  const MultilevelErrors &errors) {
    return s << "rmse=" << errors.rmse << " mse=" << errors.mse
             << " input=" << errors.input << " disc^2=" << errors.disc2;
  }
};

class EstimatorMap {
private:
  std::map<int, SingleLevelEstimator> estimators{};

  MLEstimatorConfig mlEstmConf;

  MultilevelErrors errors;

  Exponents exponents;

public:
  EstimatorMap() : mlEstmConf(MLEstimatorConfig()) {};

  explicit EstimatorMap(MLEstimatorConfig mlEstmConf)
      : mlEstmConf(std::move(mlEstmConf)), estimators(Fill(mlEstmConf)) {};

  explicit EstimatorMap(const std::map<int, WelfordAggregate> &aggregateMap,
                        MLEstimatorConfig mlEstmConf = MLEstimatorConfig())
      : mlEstmConf(std::move(mlEstmConf)), estimators(Fill(aggregateMap)) {};

  static std::map<int, SingleLevelEstimator> Fill(
      const std::map<int, WelfordAggregate> &aggregateMap);

  static std::map<int, SingleLevelEstimator> Fill(MLEstimatorConfig mlEstmConf);

  bool OpenSamples();

  double Cost() const;

  double Value() const;

#ifdef AGGREGATE_FOR_SOLUTION

  Vector MeanField() const;

  Vector SVarField() const;

#endif

  void PostSmoothing();

  void ResetSampleAmount();

  void UpdateExponentsAndErrors();

  double CostPrediction() const;

  double MemoryPrediction() const;

  double SynchronizationEfficiency() const;

  double CostWithoutSyncLosses() const;

  double SampleCounterFactor(double epsilon, double theta);

  void UpdateSampleCounterOnMap(double epsilon, double theta);

  int OptimalMOnLevel(int level, double epsilon, double theta);

  int OptimalMOnLevel(double sVarY, double meanC, double epsilon, double theta);

  bool UpdateNewLevel(double epsilon, double theta, int newLevel);

  const MultilevelErrors &GetErrors() const { return errors; };

  const Exponents &GetExponents() const { return exponents; };

  std::vector<double> GetCostPerSampleVector() const;

  std::vector<double> GetMeanQVector() const;

  std::vector<double> GetsVarQVector() const;

  std::vector<double> GetSkewQVector() const;

  std::vector<double> GetKurtQVector() const;

  std::vector<double> GetMeanYVector() const;

  std::vector<double> GetsVarYVector() const;

  std::vector<double> GetSkewYVector() const;

  std::vector<double> GetKurtYVector() const;

  std::vector<double> GetCostVector() const;

  std::vector<int> GetCommSplitVector() const;

  std::vector<int> GetLevelVector() const;

  std::vector<int> GetdMVector() const;

  std::vector<int> GetMVector() const;

  auto clear() { estimators.clear(); }

  auto end() { return estimators.end(); }

  auto rend() { return estimators.rend(); }

  auto begin() { return estimators.begin(); }

  auto rbegin() { return estimators.rbegin(); }

  auto size() const { return estimators.size(); }

  auto end() const { return estimators.end(); }

  auto rend() const { return estimators.end(); }

  auto begin() const { return estimators.begin(); }

  auto rbegin() const { return estimators.rbegin(); }

  auto operator[](int level) { return estimators[level]; }

  auto find(int level) { return estimators.find(level); }

  auto at(int level) const { return estimators.at(level); }

  auto erase(int level) { return estimators.erase(level); }

  auto insert(std::pair<int, SingleLevelEstimator> &&pair) {
    estimators.insert(pair);
  }

  friend std::ostream &operator<<(std::ostream &s, const EstimatorMap &estMap);
};

#endif  // LEVELMAPS_HPP
