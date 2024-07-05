#include "MultiLevelEstimator.hpp"
#include "Plotting.hpp"
#include <ranges>

void MultiLevelEstimator::Method() {
  if (epsilon == 0.0 && mlEstmConf.timeBudget == 0.0) MLE();
  else if (epsilon != 0.0 && mlEstmConf.timeBudget == 0.0) AMLE();
  else if (epsilon == 0.0 && mlEstmConf.timeBudget != 0.0) BMLE();
  else Exit("Configuration error")
}

void MultiLevelEstimator::MLE() {
  mout.StartBlock("MLE");
  vout(1) << "dM0=" << vec2str(estMap.GetdMVector()) << endl;

  if (!levelLoop()) Warning("No initial samples");

  mout.EndBlock(verbose == 0);
}

void MultiLevelEstimator::AMLE() {
  mout.StartBlock("AMLE");
  vout(1) << "TargetEpsilon=" << epsilon
          << " and dM0=" << vec2str(estMap.GetdMVector()) << endl;

  if (!levelLoop()) Warning("No initial samples");

  while (estMap.GetErrors().rmse >= epsilon)
    if (appendLevelIfNeeded() || updateSamplesIfNeeded())
      if (!levelLoop())
        continue;

  mout.EndBlock(verbose == 0);
}

void MultiLevelEstimator::BMLE() {
  mout.StartBlock("BMLE");
  vout(1) << "TimeBudget=" << mlEstmConf.timeBudget
          << " and dM0=" << vec2str(estMap.GetdMVector()) << endl;

  if (!levelLoop()) Warning("No initial samples")

  epsilon = mlEstmConf.eta * estMap.GetErrors().rmse;
  contData.PushBack(estMap, epsilon);

  while (!timeBudgetIsExhausted()) {
    if (appendLevelIfNeeded() || updateSamplesIfNeeded()) {
      if (expectedToFinish()) {
        if (levelLoop()) {
          contData.PushBack(estMap, epsilon);
        }
      } else if (tryUpdateEpsilon())
        continue;
      else
        break;
    }
    epsilon = mlEstmConf.eta * epsilon;
  }

  mout.EndBlock(verbose == 0);
}

bool MultiLevelEstimator::levelLoop() {
  if (!estMap.OpenSamples())
    return false;

  mout.StartBlock("LevelLoop");
  vout(1) << "Start with expected cost " << estMap.CostPrediction() << endl;
  if (verbose >= 2)
    MemoryLogger::LogMemory("Before dM=" + vec2str(estMap.GetdMVector()));
  std::map<int, std::pair<double, double>> timeMap{};
  for (auto &iter: std::ranges::reverse_view(estMap)) {
    if (iter.second.SamplesToCompute() > 0) {
      timeMap[iter.first].first = MPI_Wtime();
      iter.second.Method();
      timeMap[iter.first].second = MPI_Wtime();
    }
  }
  if (verbose >= 2)
    MemoryLogger::LogMemory("After dM=" + vec2str(estMap.GetdMVector()));

  updateData(timeMap);

  vout(2) << estMap;
  vout(2) << endl;
  vout(2) << estMap.GetExponents();
  vout(2) << endl;
  vout(1) << estMap.GetErrors();
  vout(1) << endl;
  vout(1) << "E[Q]=" << estMap.Value();
  vout(1) << " C_total=" << estMap.Cost();
  vout(1) << endl;

#ifdef AGGREGATE_FOR_SOLUTION
    if (plotting >= 1) {
      mpp::PlotData("MeanField", "MeanField." + vec2str(estMap.GetMVector()), estMap.MeanField());
      mpp::PlotData("SVarField", "SVarField." + vec2str(estMap.GetMVector()), estMap.SVarField());
    }
#endif

  mout.EndBlock(verbose == 0);

  return true;
}

void MultiLevelEstimator::updateData(
    std::map<int, std::pair<double, double>> timeMap) {
  for (auto &iter: std::ranges::reverse_view(estMap)) {
    if (iter.second.SamplesToCompute() > 0) {
      if (mlEstmConf.slEstmConf.delayParallelUpdate) {
        iter.second.UpdateParallel();
      }
      if (mlEstmConf.slEstmConf.delayTotalUpdate) {
        iter.second.UpdateTotal();
        iter.second.UpdateErrors();

      }
      if (mlEstmConf.slEstmConf.delayCostUpdate) {
        iter.second.UpdateCost(timeMap[iter.first].second -
                               timeMap[iter.first].second);
      }
    }
  }
  estMap.UpdateExponentsAndErrors();
}

void MultiLevelEstimator::ExponentResults() const {
  mout.PrintInfo("Exponents", 1,
                 PrintInfoEntry("Final alpha", estMap.GetExponents().alpha),
                 PrintInfoEntry("Final beta", estMap.GetExponents().beta),
                 PrintInfoEntry("Final gamma", estMap.GetExponents().gamma));
}

void MultiLevelEstimator::ContinuationResults() const {
  if (contData.IsEmpty()) return;
  contData.PrintInfo();
}

bool MultiLevelEstimator::timeBudgetIsExhausted() {
  if (0.95 * mlEstmConf.timeBudget <= estMap.Cost()) {
    vout(1) << "Budget is exhausted over 95%" << endl;
    return true;
  }
  return false;
}

bool MultiLevelEstimator::expectedToFinish() {
  double costPrediction = estMap.CostPrediction();
  double leftOverBudget = mlEstmConf.timeBudget - estMap.Cost();

  if (costPrediction >= leftOverBudget) {
    vout(1) << "CostPrediction=" << costPrediction
            << " leftOverBudget=" << leftOverBudget
            << " epsilon=" << epsilon << endl;
    return false;
  } else {
    attemptCtr = 0;
    if (costPrediction != 0.0)
      vout(1) << "Estimation round i=" << contData.epsilons.size()
              << " to reach epsilon=" << epsilon << endl;
    return true;
  }
}

bool MultiLevelEstimator::expectedToFitOnMemory() {
  double memPrediction = estMap.MemoryPrediction();

  if (memPrediction >= mlEstmConf.memoryBudget) {
    vout(1) << "MemoryPrediction=" << memPrediction
            << " memoryBudget=" << mlEstmConf.memoryBudget
            << " epsilon=" << epsilon << endl;
    return false;
  } else {
    return true;
  }
}

bool MultiLevelEstimator::tryUpdateEpsilon() {
  if (attemptCtr <= 30) {
    attemptCtr++;
    vout(2) << "Epsilons=" << vec2str(contData.epsilons) << " " << epsilon
            << endl;
    epsilon = 0.5 * (contData.epsilons.back() + epsilon);
    estMap.ResetSampleAmount();
    vout(2) << "Increased epsilon=" << epsilon << endl;
    return true;
  } else {
    vout(1) << "Can't increase epsilon more " << epsilon << endl;
    return false;
  }
}

bool MultiLevelEstimator::appendLevelIfNeeded() {
  if (estMap.GetErrors().disc < sqrt(1 - mlEstmConf.theta) * epsilon)
    return false;

  if (estMap.UpdateNewLevel(epsilon, mlEstmConf.theta,
                            estMap.rbegin()->first + 1)) {
    vout(1) << "Appended level=" << estMap.rbegin()->first
            << " with dM=" << estMap.rbegin()->second.SamplesToCompute()
            << endl;

    return true;
  } else {
    tryUpdateEpsilon();
    return false;
  }
}

bool MultiLevelEstimator::updateSamplesIfNeeded() {
  if (estMap.GetErrors().input < mlEstmConf.theta * epsilon * epsilon)
    return false;

  estMap.UpdateSampleCounterOnMap(epsilon, mlEstmConf.theta);

  vout(1) << "dM=" << vec2str(estMap.GetdMVector()) << endl;

  return true;
}

void MultiLevelEstimator::MultilevelResults() const {
  mout.PrintInfo("Multilevel Results", 1,
                 PrintInfoEntry("E(Qf-Qc)", estMap.GetMeanYVector()),
                 PrintInfoEntry("E(Qf)", estMap.GetMeanQVector()),
                 PrintInfoEntry("V(Qf-Qc)", estMap.GetsVarYVector()),
                 PrintInfoEntry("V(Qf)", estMap.GetsVarQVector()),
                 PrintInfoEntry("S(Qf-Qc)", estMap.GetSkewYVector()),
                 PrintInfoEntry("S(Qf)", estMap.GetSkewQVector()),
                 PrintInfoEntry("K(Qf-Qc)", estMap.GetKurtYVector()),
                 PrintInfoEntry("K(Qf)", estMap.GetKurtQVector()),
                 PrintInfoEntry("E(cost)", estMap.GetCostPerSampleVector()),
                 PrintInfoEntry("Cost on Level", estMap.GetCostVector()),
                 PrintInfoEntry("Used Levels", estMap.GetLevelVector()),
                 PrintInfoEntry("Used Samples", estMap.GetMVector()));
}

void ContinuationData::PrintInfo() const {
  mout.PrintInfo("Continuation", 1, PrintInfoEntry("Costs", costs),
                 PrintInfoEntry("Values", values),
                 PrintInfoEntry("Memory", memory),
                 PrintInfoEntry("Epsilons", epsilons),
                 PrintInfoEntry("MS Errors", msErrors),
                 PrintInfoEntry("RMS Errors", rmsErrors),
                 PrintInfoEntry("Disc Errors", discErrors),
                 PrintInfoEntry("Disc2 Errors", disc2Errors),
                 PrintInfoEntry("Input Errors", inputErrors),
                 PrintInfoEntry("CommSplits", commSplits));
}
