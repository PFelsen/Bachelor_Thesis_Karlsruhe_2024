#include "EstimatorMap.hpp"


void Exponents::ComputeAlpha(const EstimatorMap &mcMap) {
  auto _levels = mcMap.GetLevelVector();
  auto _avgsY = mcMap.GetMeanYVector();

  vector<double> levels, avgsY;
  for (int i = 1; i < _levels.size(); i++) {
    if (_avgsY[i] == 0.0) continue;
    levels.push_back(double(_levels[i]));
    avgsY.push_back(-log2(std::abs(_avgsY[i])));
  }

  alpha = ft::linearFit(levels, avgsY).first;
}

void Exponents::ComputeBeta(const EstimatorMap &mcMap) {
  auto _levels = mcMap.GetLevelVector();
  auto _varsY = mcMap.GetsVarYVector();

  vector<double> levels, varsY;
  for (int i = 1; i < _levels.size(); i++) {
    if (_varsY[i] == 0.0) continue;
    levels.push_back((double) _levels[i]);
    varsY.push_back(-log2(_varsY[i]));
  }

  beta = ft::linearFit(levels, varsY).first;
}

void Exponents::ComputeGamma(const EstimatorMap &mcMap) {
  auto _levels = mcMap.GetLevelVector();
  auto _avgsCost = mcMap.GetCostPerSampleVector();

  vector<double> levels, avgsCost;
  for (int i = 1; i < _levels.size(); i++) {
    if (_avgsCost[i] == 0.0) continue;
    levels.push_back((double) _levels[i]);
    avgsCost.push_back(log2(_avgsCost[i]));
  }

  gamma = ft::linearFit(levels, avgsCost).first;
}

void WelfordAggregate::UpdateErrors() {
  errors.EstimateErrors(*this);
}


void MultilevelErrors::EstimateErrors(const EstimatorMap &mcMap) {
  disc = EstimateNumeric(mcMap);
  disc2 = disc * disc;
  input = EstimateStochastic(mcMap);
  mse = input + disc2;
  rmse = sqrt(mse);
}

double MultilevelErrors::EstimateNumeric(const EstimatorMap &mcMap) {
  auto _levels = mcMap.GetLevelVector();
  auto _means = mcMap.GetMeanYVector();
  double alpha = mcMap.GetExponents().alpha;

  double err = 0.0;
  int integer = (int) mcMap.size() - 2;
  for (int i = 1; i < _levels.size(); i++) {
    if (_means[i] == 0.0) continue;
    err = max(std::abs(_means[i]) / (pow(2.0, alpha) - 1) / (pow(2.0, integer)), err);
    integer--;
  }

  return err;
}

double MultilevelErrors::EstimateStochastic(const EstimatorMap &mcMap) {
  auto _ctrs = mcMap.GetMVector();
  auto _sVars = mcMap.GetsVarYVector();

  double err = 0.0;
  for (int i = 0; i < _ctrs.size(); i++) {
    if (_ctrs[i] == 0.0) break;
    err += _sVars[i] / double(_ctrs[i]);
  }
  return err;
}

void EstimatorMap::UpdateExponentsAndErrors() {
  exponents.ComputeAlpha(*this);
  exponents.ComputeBeta(*this);
  exponents.ComputeGamma(*this);
  errors.EstimateErrors(*this);
}

void EstimatorMap::ResetSampleAmount() {
  for (auto &[level, estimator]: estimators)
    estimator.UpdateSampleAmount(0);
}

double EstimatorMap::SampleCounterFactor(double epsilon, double theta) {
  double factor = 0.0;
  for (auto &[level, estimator]: estimators)
    factor += sqrt(estimator.GetAggregate().Y.GetSVar() * estimator.GetAggregate().costPerSample);
  factor *= pow(sqrt(theta) * epsilon, -2);
  return factor;
}

int EstimatorMap::OptimalMOnLevel(int level, double epsilon, double theta) {
  double sVarY = estimators[level].GetAggregate().Y.GetSVar();
  double meanC = estimators[level].GetAggregate().costPerSample;

  return OptimalMOnLevel(sVarY, meanC, epsilon, theta);
}

int EstimatorMap::OptimalMOnLevel(double sVarY, double meanC, double epsilon, double theta) {
  return (int) (ceil(SampleCounterFactor(epsilon, theta) * sqrt(sVarY / meanC)));
}

void EstimatorMap::UpdateSampleCounterOnMap(double epsilon, double theta) {
  for (auto &[level, estimator]: estimators)
    estimator.UpdateSampleAmount(
        OptimalMOnLevel(level, epsilon, theta) - estimator.GetAggregate().ctr.M
    );
  PostSmoothing();
}

void EstimatorMap::PostSmoothing() {
  for (auto &[level, estimator]: estimators) {
    if (level == estimators.begin()->first) continue;
    if (level == estimators.rbegin()->first) break;
    int expectedMonLevel = estimator.GetAggregate().ctr.M + estimator.SamplesToCompute();
    int expectedMafterLevel = estimators[level + 1].GetAggregate().ctr.M
                              + estimators[level + 1].SamplesToCompute();
    int expectedMbeforeLevel = estimators[level - 1].GetAggregate().ctr.M
                               + estimators[level - 1].SamplesToCompute();
    if (expectedMonLevel <= expectedMafterLevel)
      estimator.UpdateSampleAmount(
          ((int) pow(2, ((log2(expectedMafterLevel) + log2(expectedMbeforeLevel)) / 2.0)))
          - estimator.GetAggregate().ctr.M);
  }
}

bool EstimatorMap::UpdateNewLevel(double epsilon, double theta, int newLevel) {
  double costGrowth = pow(2.0, exponents.gamma);
  double varianceReduction = pow(2.0, -exponents.beta);
  double sVarY = varianceReduction * estimators[newLevel - 1].GetAggregate().Y.GetSVar();
  double costPerSample = costGrowth * estimators[newLevel - 1].GetAggregate().costPerSample;

  int optimalMOnLevel = OptimalMOnLevel(sVarY, costPerSample, epsilon, theta);
  double costPrediction = costPerSample * optimalMOnLevel;
  double leftOverBudget = mlEstmConf.timeBudget - Cost();

  if (costPrediction >= leftOverBudget) return false;

  this->insert({newLevel, SingleLevelEstimator(mlEstmConf.slEstmConf.WithInitLevel(newLevel).
      WithOnlyFine(false).WithInitSamples(optimalMOnLevel)
  )});

  estimators[newLevel].UpdateSVarYAndCostPerSample(sVarY, costPerSample); // Todo is this one needed?
  return true;
}

std::map<int, SingleLevelEstimator> EstimatorMap::Fill(MLEstimatorConfig _mlEstmConf) {
  std::map<int, SingleLevelEstimator> _estimators{};
  for (auto level_index = 0; level_index < _mlEstmConf.initLevels.size(); level_index++) {
    _estimators.insert({_mlEstmConf.initLevels[level_index], (
        SingleLevelEstimator(
            _mlEstmConf.slEstmConf.WithInitSamples(_mlEstmConf.initSamples[level_index]).
                WithInitLevel(_mlEstmConf.initLevels[level_index]).WithEpsilon(0.0).
                WithOnlyFine((level_index == 0))))});
  }
  return _estimators;
}

std::map<int, SingleLevelEstimator>
EstimatorMap::Fill(const std::map<int, WelfordAggregate> &aggregateMap) {
  std::map<int, SingleLevelEstimator> _estimators{};
  for (auto const &[level, aggregate]: aggregateMap) {
    _estimators.insert({level, SingleLevelEstimator(
        SLEstimatorConfig().WithInitLevel(level), aggregate)});
  }
  return _estimators;
}

double EstimatorMap::Cost() const {
  double cost = 0.0;
  for (auto &[level, estimator]: estimators)
    cost += estimator.GetAggregate().cost;
  return cost;
}

double EstimatorMap::Value() const {
  double value = 0.0;
  for (auto &[level, estimator]: estimators)
    value += estimator.GetAggregate().Y.GetMean();
  return value;
}

#ifdef AGGREGATE_FOR_SOLUTION

#include "Transfers.hpp"

Vector EstimatorMap::MeanField() const {
  auto coarseMeanField = std::make_shared<Vector>(
      estimators.at(mlEstmConf.initLevels[0]).GetAggregate().V.GetMean()
  );
  const auto levels = GetLevelVector();
  for (int i = 1; i < estimators.size(); i++) {
    Vector transferCoarseField(0.0, estimators.at(levels[i]).GetAggregate().V.GetMean());
    Vector fineField(estimators.at(levels[i]).GetAggregate().V.GetMean());
    auto transfer = GetTransfer(*coarseMeanField, transferCoarseField, "Matrix");
    transfer->Prolongate(*coarseMeanField, transferCoarseField);
    transferCoarseField += fineField;
    coarseMeanField = std::make_shared<Vector>(transferCoarseField);
  }
  return *coarseMeanField;
}

Vector EstimatorMap::SVarField() const {
  const auto levels = GetLevelVector();
  const auto ctrs = GetMVector();
  auto coarseSVarField = std::make_shared<Vector>(
      estimators.at(mlEstmConf.initLevels[0]).GetAggregate().V.GetSVar()
  );
  *coarseSVarField *= 1.0 / ctrs[0];
  for (int i = 1; i < estimators.size(); i++) {
    Vector transferCoarseField(0.0, estimators.at(levels[i]).GetAggregate().V.GetSVar());
    Vector fineField(estimators.at(levels[i]).GetAggregate().V.GetSVar());
    auto transfer = GetTransfer(*coarseSVarField, transferCoarseField, "Matrix");
    transfer->Prolongate(*coarseSVarField, transferCoarseField);
    fineField *= 1.0 / ctrs[i];
    transferCoarseField += fineField;
    coarseSVarField = std::make_shared<Vector>(transferCoarseField);
  }
  return *coarseSVarField;
}

#endif

double EstimatorMap::CostPrediction() const {
  double costPrediction = 0.0;
  for (auto const &[level, estimator]: estimators)
    costPrediction += estimator.SamplesToCompute() * estimator.GetAggregate().costPerSample;
  return costPrediction;
}

double EstimatorMap::MemoryPrediction() const {
  double memoryConstant = 1.0;
  double memoryPrediction = 0.0;
  for (auto const &[level, estimator]: estimators)
    memoryPrediction += memoryConstant * estimator.CellsOn(LevelPair{level, -1, 0, estimator.CommSplit()});
  return memoryPrediction;
}



double EstimatorMap::CostWithoutSyncLosses() const {
  double costWithoutSyncLosses = 0.0;
  for (auto &[level, estimator]: estimators)
    costWithoutSyncLosses += estimator.GetAggregate().C.GetMean() * estimator.GetAggregate().ctr.M;
  return costWithoutSyncLosses;
}

double EstimatorMap::SynchronizationEfficiency() const {
  return CostWithoutSyncLosses() / Cost();
}

bool EstimatorMap::OpenSamples() {
  for (auto const &[level, estimator]: estimators)
    if (estimator.SamplesToCompute() > 0)
      return true;
  return false;
}

std::ostream &operator<<(std::ostream &s, const EstimatorMap &estMap) {
  return s << "Level=" << vec2str(estMap.GetLevelVector()) << endl
           << "M_lvl=" << vec2str(estMap.GetMVector()) << endl
           << "MeanY=" << vec2str(estMap.GetMeanYVector()) << endl
           << "MeanQ=" << vec2str(estMap.GetMeanQVector()) << endl
           << "SVarY=" << vec2str(estMap.GetsVarYVector()) << endl
           << "SVarQ=" << vec2str(estMap.GetsVarQVector()) << endl
           << "SkewY=" << vec2str(estMap.GetSkewYVector()) << endl
           << "KurtY=" << vec2str(estMap.GetKurtYVector()) << endl
           << "MeanC=" << vec2str(estMap.GetCostPerSampleVector()) << endl
           << "TCost=" << vec2str(estMap.GetCostPerSampleVector()) << endl;
}

std::ostream &operator<<(std::ostream &s, const MLEstimatorConfig &conf) {
  return s << DOUT(conf.eta)
           << DOUT(conf.theta)
           << DOUT(conf.epsilon)
           << DOUT(conf.timeBudget)
           << "conf.initLevels:" << vec2str(conf.initLevels)
           << " conf.initSamples:" << vec2str(conf.initSamples);
}

std::vector<int> EstimatorMap::GetdMVector() const {
  std::vector<int> dMVector;
  for (auto const &[level, estimator]: estimators)
    dMVector.push_back(estimator.SamplesToCompute());
  return dMVector;
}

std::vector<int> EstimatorMap::GetCommSplitVector() const {
  std::vector<int> minCommSplitVector;
  for (auto const &[level, estimator]: estimators)
    minCommSplitVector.push_back(estimator.CommSplit());
  return minCommSplitVector;
}

std::vector<int> EstimatorMap::GetLevelVector() const {
  std::vector<int> levelVector;
  for (auto const &[level, estimator]: estimators)
    levelVector.push_back(estimator.Level());
  return levelVector;
}

std::vector<double> EstimatorMap::GetMeanQVector() const {
  std::vector<double> qVector;
  for (auto const &[level, estimator]: estimators)
    qVector.push_back(estimator.GetAggregate().Q.GetMean());
  return qVector;
}

std::vector<double> EstimatorMap::GetsVarQVector() const {
  std::vector<double> qVector;
  for (auto const &[level, estimator]: estimators)
    qVector.push_back(estimator.GetAggregate().Q.GetSVar());
  return qVector;
}

std::vector<double> EstimatorMap::GetSkewQVector() const {
  std::vector<double> qVector;
  for (auto const &[level, estimator]: estimators)
    qVector.push_back(estimator.GetAggregate().Q.GetSkew());
  return qVector;
}

std::vector<double> EstimatorMap::GetKurtQVector() const {
  std::vector<double> qVector;
  for (auto const &[level, estimator]: estimators)
    qVector.push_back(estimator.GetAggregate().Q.GetKurt());
  return qVector;
}

std::vector<double> EstimatorMap::GetMeanYVector() const {
  std::vector<double> yVector;
  for (auto const &[level, estimator]: estimators)
    yVector.push_back(estimator.GetAggregate().Y.GetMean());
  return yVector;
}

std::vector<double> EstimatorMap::GetsVarYVector() const {
  std::vector<double> yVector;
  for (auto const &[level, estimator]: estimators)
    yVector.push_back(estimator.GetAggregate().Y.GetSVar());
  return yVector;
}

std::vector<double> EstimatorMap::GetSkewYVector() const {
  std::vector<double> yVector;
  for (auto const &[level, estimator]: estimators)
    yVector.push_back(estimator.GetAggregate().Y.GetSkew());
  return yVector;
}

std::vector<double> EstimatorMap::GetKurtYVector() const {
  std::vector<double> yVector;
  for (auto const &[level, estimator]: estimators)
    yVector.push_back(estimator.GetAggregate().Y.GetKurt());
  return yVector;
}

std::vector<double> EstimatorMap::GetCostPerSampleVector() const {
  std::vector<double> costVector;
  for (auto const &[level, estimator]: estimators)
    costVector.push_back(estimator.GetAggregate().costPerSample);
  return costVector;
}

std::vector<double> EstimatorMap::GetCostVector() const {
  std::vector<double> costVector;
  for (auto const &[level, estimator]: estimators)
    costVector.push_back(estimator.GetAggregate().cost);
  return costVector;
}

std::vector<int> EstimatorMap::GetMVector() const {
  std::vector<int> mVector;
  for (auto const &[level, estimator]: estimators)
    mVector.push_back(estimator.GetAggregate().ctr.M);
  return mVector;
}
