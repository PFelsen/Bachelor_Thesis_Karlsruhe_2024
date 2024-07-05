#ifndef MPP_MAIN_HPP
#define MPP_MAIN_HPP

#include "SingleExecution.hpp"
#include "ConvergenceStudy.hpp"
#include "AdaptiveStudy.hpp"


SingleExecution SingleExperiment() {
  auto singleExecution = SingleExecution(PDESolverConfig());
  singleExecution.Method();
  return singleExecution;
}

std::shared_ptr<ConvergenceStudy> ConvergenceExperiment() {
  auto convergenceStudy = CreateConvergenceStudy(PDESolverConfig());
  convergenceStudy->Method();
  return convergenceStudy;
}

std::shared_ptr<AdaptiveStudy> AdaptiveExperiment() {
  auto adaptiveStudy = CreateAdaptiveStudy(PDESolverConfig(), AdaptiveConfig());
  adaptiveStudy->Method();
  return adaptiveStudy;
}

#if USE_SPACETIME

#include "FWISTExperiment.hpp"

FWISTExperiment FWIExperimentST() {
  auto fwiExperiment = FWISTExperiment(PDESolverConfig());
  fwiExperiment.Method();
  return fwiExperiment;
}

#endif

#if BUILD_UQ

#include "StochasticGradientDescent.hpp"
#include "MultiLevelEstimator.hpp"
#include "SingleLevelEstimator.hpp"
#include "SequentialMonteCarlo.hpp"
#include "Random.hpp"

StochasticGradientDescent SGDExperiment() {
  Random::Initialize();
  auto sgd = StochasticGradientDescent();
  sgd.Method();
  return sgd;
}

SequentialMonteCarlo SMCExperiment() {
  Random::Initialize();
  auto estimator = SequentialMonteCarlo(SequentialMonteCarloConfig());
  estimator.Method();
  estimator.EstimatorResults();
  return estimator;
}

SingleLevelEstimator MCExperiment() {
  Random::Initialize();
  auto estimator = SingleLevelEstimator(SLEstimatorConfig());
  estimator.Method();
  estimator.EstimatorResults();
  return estimator;
}

MultiLevelEstimator MLMCExperiment() {
  MemoryLogger::Tare();
  mout.StartBlock("MLMCExperiment");
  mout << "Start" << endl;
  Random::Initialize();
  auto estimator = MultiLevelEstimator(MLEstimatorConfig());
  estimator.Method();
  mout.EndBlock();
  mout << endl;
  estimator.MultilevelResults();
  estimator.ExponentResults();
  estimator.EstimatorResults();
  estimator.ContinuationResults();
  MemoryLogger::PrintInfo();
  return estimator;
}

// Todo: Setup SC Experiment.

#endif

#endif //MPP_MAIN_HPP
