#include "SingleLevelEstimator.hpp"
#include "MemoryLogger.hpp"
#include "Plotting.hpp"
#include <utility>

void SingleLevelEstimator::Method() {
  mout.StartBlock(std::string(slEstmConf.parallel ? "Parallel" : "Seriell")
                  + "Estimator" + to_string(initLevel));

  double start = MPI_Wtime();
  if (slEstmConf.epsilon == 0.0) {
    vout(1) << "dM=" << dM << " dMcomm=" << dMcomm << " commSplit="
            << commSplit << " " << aggregate.ctr << endl;
  } else {
    vout(1) << "epsilon=" << slEstmConf.epsilon << endl;
  }

  MemoryLogger::LogMemory();

  while (dM > 0) {
    while (dMcomm > 0) {
      auto fineSolution = msFEM->RunProtocol(SampleID(initLevel, Index(), false, commSplit));
      auto coarseSolution = slEstmConf.onlyFine ?
                            msFEM->EmptyProtocol(SampleID(initLevel, Index(), false, commSplit)) :
                            msFEM->RunProtocol(SampleID(initLevel, Index(), true, commSplit));
      UpdateOnComm(fineSolution, coarseSolution);
      dMcomm--;
    }
    MemoryLogger::LogMemory();
    if (!slEstmConf.delayParallelUpdate) UpdateParallel();
    MemoryLogger::LogMemory();
    if (!slEstmConf.delayTotalUpdate) UpdateTotal();
    MemoryLogger::LogMemory();
    if (!slEstmConf.delayTotalUpdate) UpdateErrors();
    updateSamplesIfNeeded();
  }
  MemoryLogger::LogMemory();

  if (!slEstmConf.delayCostUpdate)
    UpdateCost(MPI_Wtime() - start);

  mout.EndBlock(verbose == 0);
}

void SingleLevelEstimator::updateSamplesIfNeeded() {
  if (slEstmConf.epsilon == 0.0) {
    dM = 0;
  } else {
    if (aggregate.errors.disc < slEstmConf.epsilon) {
      if (aggregate.errors.rmse > slEstmConf.epsilon) {
        vout(2) << "UpdateOnComm sample amount" << endl;
        int optimalM = (int) (ceil(aggregate.Q.GetSVar() / (pow(slEstmConf.epsilon, 2))));
        UpdateSampleAmount(optimalM - aggregate.ctr.M);
      } else {
        dM = 0;
      }
    } else {
      Warning("Numerical error is too big to reach error bound")
    }
  }
}

void SingleLevelEstimator::UpdateSampleAmount(int requestedSamples) {
  if (requestedSamples <= 0) {
    this->commSplit = 0;
    this->dMcomm = 0;
    this->dM = 0;
  } else {
    // Variance is otherwise not computable
    if (requestedSamples == 1)
      requestedSamples++;
    if (slEstmConf.parallel) {
      for (int i = 0; i < ceil(log2(PPM->Size(0))) + 1; i++) {
        if (requestedSamples >= PPM->Size(i)) {
          this->dM = requestedSamples;
          int mod = requestedSamples % PPM->Size(i);
          if (mod != 0)
            this->dM += (PPM->Size(i) - mod);
          this->dMcomm = ceil((double) requestedSamples / PPM->Size(i));
          this->commSplit = PPM->MaxCommSplit() - i;
          break;
        }
      }
    } else {
      this->commSplit = 0;
      this->dMcomm = requestedSamples;
      this->dM = requestedSamples;
    }
  }
}

int SingleLevelEstimator::Index() const {
  int index =
      aggregate.ctr.M + aggregate.ctr.Mcomm +
      (dM / PPM->NumberOfCommunicators(commSplit)) * PPM->Color(commSplit);
  return index;
}

void SingleLevelEstimator::UpdateOnComm(const SampleSolution &fSolution,
                                        const SampleSolution &cSolution) {
  aggregate.UpdateOnComm(fSolution, cSolution);
  if (verbose > 2) {
    mout << endl << "Data after UpdateOnComm()" << endl;
    pout << "Q: " << aggregate.Q.GetComm() << "Y: " << aggregate.Y.GetComm()
         << "C: " << aggregate.C.GetComm() << endl;
  }
#ifdef AGGREGATE_FOR_SOLUTION
  if (plotting >= 3) {
    std::string nameSubStr = std::to_string(aggregate.U.GetComm().mean.SpaceLevel()) +
                             ".MeanComm." + std::to_string(aggregate.ctr.Mcomm);
    mpp::PlotData("U", "U." + nameSubStr, aggregate.U.GetComm().mean);
    mpp::PlotData("V", "V." + nameSubStr, aggregate.V.GetComm().mean);
  }
#endif
}

void SingleLevelEstimator::UpdateParallel() {
  aggregate.UpdateParallel(commSplit);
  msFEM->ClearMeshesOnMeshIndex(LevelPair{initLevel, -1, 0, commSplit});
  if (!slEstmConf.onlyFine)
    msFEM->ClearMeshesOnMeshIndex(LevelPair{initLevel - 1, -1, 0, commSplit});

  if (verbose > 2) {
    mout << endl << "Data after UpdateParallel()" << endl;
    pout << "Q: " << aggregate.Q.GetPara() << "Y: " << aggregate.Y.GetPara()
         << "C: " << aggregate.C.GetPara() << endl;
  }
#ifdef AGGREGATE_FOR_SOLUTION
  if (plotting >= 2) {
    std::string nameSubStr = std::to_string(aggregate.U.GetPara().mean.SpaceLevel()) +
                             ".MeanPara." + std::to_string(aggregate.ctr.Mpara);
    mpp::PlotData("U", "U." + nameSubStr, aggregate.U.GetPara().mean);
    mpp::PlotData("V", "V." + nameSubStr, aggregate.V.GetPara().mean);
  }
#endif
}

void SingleLevelEstimator::UpdateTotal() {
  aggregate.UpdateTotal();
  if (verbose > 2) {
    mout << endl << "Data after UpdateTotal()" << endl;
    pout << "Q: " << aggregate.Q.GetTotal() << "Y:" << aggregate.Y.GetTotal()
         << "C: " << aggregate.C.GetTotal() << endl;
  }
#ifdef AGGREGATE_FOR_SOLUTION
  if (plotting >= 1) {
    std::string nameSubStr = std::to_string(aggregate.U.GetTotal().mean.SpaceLevel()) +
                             ".MeanTotal." + std::to_string(aggregate.ctr.M);
    mpp::PlotData("U", "U." + nameSubStr, aggregate.U.GetTotal().mean);
    mpp::PlotData("V", "V." + nameSubStr, aggregate.V.GetTotal().mean);
  }
#endif
}

void SingleLevelEstimator::UpdateErrors() {
  if (slEstmConf.delayTotalUpdate) return;
  aggregate.UpdateErrors();
  if (verbose > 2) {
    pout << aggregate.errors << endl;
  }
}


void SingleLevelEstimator::UpdateCost(double newRoundCost) {
  aggregate.UpdateCost(newRoundCost);
  if (verbose > 2) {
    pout << DOUT(aggregate.costPerSample) << " " << DOUT(aggregate.C.GetMean())
         << endl;
    pout << DOUT(aggregate.cost) << " "
         << DOUT(aggregate.C.GetMean() * aggregate.ctr.M) << endl;
  }
}

void SingleLevelEstimator::UpdateSVarYAndCostPerSample(
    double sVarY, double expectedCostPerSample) {
  aggregate.UpdateSVarYAndCostPerSample(sVarY, expectedCostPerSample);
}

SLEstimatorConfig SLEstimatorConfig::WithInitLevel(int level) {
  msFEMConfig.pdeSolverConf = msFEMConfig.pdeSolverConf.WithLevel(level);
  initLevel = level;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithEpsilon(double targetEpsilon) {
  epsilon = targetEpsilon;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithInitSamples(int samples) {
  initSamples = samples;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithOnlyFine(bool onlyComputeOnFine) {
  onlyFine = onlyComputeOnFine;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithParallel(bool parallelEstimator) {
  parallel = parallelEstimator;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithMultiSampleFEMConfig(const MultiSampleFEMConfig &conf) {
  msFEMConfig = conf;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithDelayTotalUpdate(bool delayTotUpdate) {
  delayTotalUpdate = delayTotUpdate;
  return *this;
}

SLEstimatorConfig SLEstimatorConfig::WithDelayParallelUpdate(bool delayParaUpdate) {
  delayParallelUpdate = delayParaUpdate;
  return *this;
}

std::ostream &operator<<(std::ostream &s, const SLEstimatorConfig &conf) {
  return s << DOUT(conf.initLevel)
           << DOUT(conf.initSamples)
           << DOUT(conf.epsilon)
           << DOUT(conf.parallel)
           << DOUT(conf.onlyFine)
           << DOUT(conf.delayCostUpdate)
           << DOUT(conf.delayTotalUpdate)
           << DOUT(conf.delayParallelUpdate) << endl;
}
